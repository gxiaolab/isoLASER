#!/usr/bin/env python
import subprocess as sp
import sys, re, copy, os, codecs, gzip
from argparse import ArgumentParser, Action as ArgParseAction
from collections import OrderedDict, Counter, defaultdict
import pysam
import numpy as np



def parse_coordinates(c):
        c = c.replace(",", "")
        chr = c.split(":")[0]
        start, end = c.split(":")[1].split("-")
        # Convert to 0-based
        start, end = int(start) - 1, int(end)
        return chr, start, end



def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):

        # Match
        if CIGAR_op == "M":
                for i in range(pos, pos + CIGAR_len):
                        if i < start or i >= end:
                                continue
                        ind = i - start
                        a[ind] += 1

        # Insertion or Soft-clip
        if CIGAR_op == "I" or CIGAR_op == "S":
                return pos

        # Deletion
        if CIGAR_op == "D":
                pass

        # Junction
        if CIGAR_op == "N":
                don = pos
                acc = pos + CIGAR_len
                if don > start and acc < end:
                        junctions[(don,acc)] = junctions.setdefault((don,acc), 0) + 1

        pos = pos + CIGAR_len

        return pos


def flip_read(s, samflag):
        if s == "NONE" or s == "SENSE":
                return 0
        if s == "ANTISENSE":
                return 1
        if s == "MATE1_SENSE":
                if int(samflag) & 64:
                        return 0
                if int(samflag) & 128:
                        return 1
        if s == "MATE2_SENSE":
                if int(samflag) & 64:
                        return 1
                if int(samflag) & 128:
                        return 0


def read_bam(f, c, s):

        chr, start, end = parse_coordinates(c)

        # Initialize coverage array and junction dict
        a = {"+" : [0] * (end - start)}
        junctions = {"+": OrderedDict()}
        if s != "NONE":
                a["-"] = [0] * (end - start)
                junctions["-"] = OrderedDict()

        samfile = pysam.AlignmentFile(f)
        txs = []
        for read in samfile.fetch(chr, start, end):

                # Move forward if read is unmapped
                if read.is_unmapped:
                    continue

                samflag, read_start, CIGAR = read.flag, read.reference_start+1, read.cigarstring

                # Ignore reads with more exotic CIGAR operators
                if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
                        continue

                read_strand = ["+", "-"][flip_read(s, samflag) ^ bool(int(samflag) & 16)]
                if s == "NONE": read_strand = "+"

                CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
                CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

                pos = read_start
                for n, CIGAR_op in enumerate(CIGAR_ops):
                        CIGAR_len = int(CIGAR_lens[n])
                        pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a[read_strand], junctions[read_strand])
               
                txs.append(read.get_tag("ZT"))

        samfile.close()
        return a, junctions, txs



def read_bam_input(f, overlay, color, label):
        with codecs.open(f, encoding='utf-8') as openf:
                for line in openf:
                        line_sp = line.strip().split("\t")
                        bam = line_sp[1]
                        overlay_level = line_sp[overlay-1] if overlay else None
                        color_level = line_sp[color-1] if color else None
                        label_text = line_sp[label-1] if label else None 
                        yield line_sp[0], bam, overlay_level, color_level, label_text


def prepare_for_R(a, junctions, c, m, f):

        _, start, _ = parse_coordinates(args.coordinates)

        # Convert the array index to genomic coordinates
        x = list(i+start for i in range(len(a)))
        y = a

        # Arrays for R
        dons, accs, yd, ya, counts = [], [], [], [], []

        # Prepare arrays for junctions (which will be the arcs)
        if junctions:
                max_n = max(junctions.values())
        else:
                max_n = np.nan
        for (don, acc), n in junctions.items():

                # Do not add junctions with less than defined coverage
                if n < m or (n/max_n < f):
                        continue

                dons.append(don)
                accs.append(acc)
                counts.append(n)

                yd.append( a[ don - start -1 ])
                ya.append( a[ acc - start +1 ])

        return x, y, dons, accs, yd, ya, counts


def read_gtf(f, c, tx_id_list):

        tx_id_dict = Counter(tx_id_list)
        n = sum(tx_id_dict.values())

        transcript_frequency = {}
        for k, v in tx_id_dict.items():
                transcript_frequency[f'"{k}"'] = float(v/n)

        exons = OrderedDict()
        transcripts = OrderedDict()
        chrom, start, end = parse_coordinates(c)
        end = end -1
        with gzip.open(f, 'rt') if f.endswith(".gz") else open(f) as openf:
                for line in openf:
                        if line.startswith("#"):
                                continue
                        el_chr, _, el, el_start, el_end, _, strand, _, tags = line.strip().split("\t")
                        if el_chr != chrom:
                                continue
                        if el not in ("transcript", "exon"):
                                continue
                        try:
                                transcript_id = re.findall('transcript_id ("[^"]+")', tags)[0]
                        except KeyError:
                                print("ERROR: 'transcript_id' attribute is missing in the GTF file.")
                                exit(1)

                        if transcript_id not in transcript_frequency.keys():
                                continue
                        
                        if transcript_frequency[transcript_id] >= 0.05 or "ENST" in transcript_id:

                                el_start, el_end = int(el_start) -1, int(el_end)
                                strand = '"' + strand + '"'
                                if el == "transcript":
                                        if (el_end > start and el_start < end):
                                                transcripts[transcript_id] = max(start, el_start), min(end, el_end), strand
                                        continue
                                if el == "exon":
                                        if (start < el_start < end or start < el_end < end):
                                                exons.setdefault(transcript_id, []).append((max(el_start, start), min(end, el_end), strand))

        return transcripts, exons


def make_introns(transcripts, exons, intersected_introns=None):
        new_transcripts = copy.deepcopy(transcripts)
        new_exons = copy.deepcopy(exons)
        introns = OrderedDict()
        if intersected_introns:
                for tx, (tx_start,tx_end,strand) in new_transcripts.items():
                        total_shift = 0
                        for a,b in intersected_introns:
                                l = b - a
                                shift = l - int(l**0.7)
                                total_shift += shift
                                for i, (exon_start,exon_end,strand) in enumerate(exons.get(tx,[])):
                                        new_exon_start, new_exon_end = new_exons[tx][i][:2]
                                        if a < exon_start:
                                                if b > exon_end:
                                                        if i ==  len(exons[tx])-1:
                                                                total_shift = total_shift - shift + (exon_start - a)*(1-int(l**-0.3))
                                                        shift = (exon_start - a)*(1-int(l**-0.3))
                                                        new_exon_end = new_exons[tx][i][1] - shift
                                                new_exon_start = new_exons[tx][i][0] - shift
                                        if b <= exon_end:
                                                new_exon_end = new_exons[tx][i][1] - shift
                                        new_exons[tx][i] = (new_exon_start,new_exon_end,strand)
                        tx_start = min(tx_start, sorted(new_exons.get(tx, [[sys.maxsize]]))[0][0])
                        new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)

        for tx, (tx_start,tx_end,strand) in new_transcripts.items():
                intron_start = tx_start
                ex_end = 0
                for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
                        intron_end = ex_start
                        if tx_start < ex_start:
                                introns.setdefault(tx, []).append((intron_start, intron_end, strand))
                        intron_start = ex_end
                if tx_end > ex_end:
                        introns.setdefault(tx, []).append((intron_start, tx_end, strand))
        d = {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}
        return d


def gtf_for_ggplot(annotation, coords, mi_file, start, end, arrow_bins):
        arrow_space = int((end - start)/arrow_bins)
        chrom, start, end = parse_coordinates(coords)

        s = """

        # data table with exons
        ann_list = list(
                "exons" = data.table(),
                "introns" = data.table()
        )
	print(getwd())

        df <- read.csv("%(mi)s", sep = '\t', header = F)
	names(df) <- c("chrom", "start", "end", "gene", "dpsi")
        df = df[df$chrom == "%(chr)s",]

        """ %({"mi" : mi_file, "chr" : chrom})

        if annotation["exons"]:

                s += """
                ann_list[['exons']] = data.table(
                        tx = rep(c(%(tx_exons)s), c(%(n_exons)s)),
                        start = c(%(exon_s)s),
                        end = c(%(exon_e)s),
                        strand = c(%(strand)s)
                )
                """ %({
                "tx_exons" : ",".join(annotation["exons"].keys()),
                "n_exons" :  ",".join(map(str, map(len, annotation["exons"].values()))),
                "exon_s" : ",".join(map(str, (v[0] for vs in annotation["exons"].values() for v in vs))),
                "exon_e" : ",".join(map(str, (v[1] for vs in annotation["exons"].values() for v in vs))),
                "strand" : ",".join(map(str, (v[2] for vs in annotation["exons"].values() for v in vs))),
                })

        if annotation["introns"]:

                s += """
                ann_list[['introns']] = data.table(
                        tx = rep(c(%(tx_introns)s), c(%(n_introns)s)),
                        start = c(%(intron_s)s),
                        end = c(%(intron_e)s),
                        strand = c(%(strand)s)
                )
                # Create data table for strand arrows
                txarrows = data.table()
                introns = ann_list[['introns']]
                # Add right-pointing arrows for plus strand
                if ("+" %%in%% introns$strand && nrow(introns[strand=="+" & end-start>500, ]) > 0) {
                        txarrows = rbind(
                                txarrows,
                                introns[strand=="+" & end-start>500, list(
                                        seq(start+4,end,by=%(arrow_space)s)-1,
                                        seq(start+4,end,by=%(arrow_space)s)
                                        ), by=.(tx,start,end)
                                ]
                        )
                }
                # Add left-pointing arrows for minus strand
                if ("-" %%in%% introns$strand && nrow(introns[strand=="-" & end-start>500, ]) > 0) {
                        txarrows = rbind(
                                txarrows,
                                introns[strand=="-" & end-start>500, list(
                                        seq(start,max(start+1, end-4), by=%(arrow_space)s),
                                        seq(start,max(start+1, end-4), by=%(arrow_space)s)-1
                                        ), by=.(tx,start,end)
                                ]
                        )
                }
                """ %({
                        "tx_introns": ",".join(annotation["introns"].keys()),
                        "n_introns": ",".join(map(str, map(len, annotation["introns"].values()))),
                        "intron_s" : ",".join(map(str, (v[0] for vs in annotation["introns"].values() for v in vs))),
                        "intron_e" : ",".join(map(str, (v[1] for vs in annotation["introns"].values() for v in vs))),
                        "strand"   : ",".join(map(str, (v[2] for vs in annotation["introns"].values() for v in vs))),
                        "arrow_space" : 5000,
                })

        s += """

        gtfp = ggplot()
        if (length(ann_list[['introns']]) > 0) {
                gtfp = gtfp + geom_segment(data=ann_list[['introns']], aes(x=start, xend=end, y=tx, yend=tx), linewidth=0.3)
                gtfp = gtfp + geom_segment(data=txarrows, aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")), linewidth = 0.3)
        }
        if (length(ann_list[['exons']]) > 0) {
	        gtfp = gtfp + geom_rect(data=df, fill='orange', alpha = 0.3, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf)
                gtfp = gtfp + geom_segment(data=ann_list[['exons']], aes(x=start, xend=end, y=tx, yend=tx), linewidth=5, alpha=1)

        }
        gtfp = gtfp + scale_y_discrete(expand=c(0,0.5))
        gtfp = gtfp + scale_x_continuous(expand=c(0,0.25))
        gtfp = gtfp + coord_cartesian(xlim = c(%s,%s))
        gtfp = gtfp + labs(y=NULL)
        gtfp = gtfp + theme_bw() #+ theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
        """ %(start, end)

        return s


def setup_R_script(label_dict):
        s = """
        library(ggplot2)
        library(grid)
        library(gridExtra)
        library(data.table)
        library(gtable)
        library(patchwork)
        library(RColorBrewer)

        scale_lwd = function(r) {
                lmin = 0.1
                lmax = 4
                return( r*(lmax-lmin)+lmin )
        }
     
        theme_update(
                plot.margin = unit(c(15,15,15,15), "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(linewidth=0.5),
                axis.title.x = element_blank(),
                axis.title.y = element_text(angle=0, vjust=0.5)
        )

        labels = list(%(labels)s)

        density_list = list()
        junction_list = list()

        """ %({ 'labels': ",".join(('"%s"="%s"' %(idx,lab) for idx,lab in label_dict.items())),
        })
        return s

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.

def mean(lst):
        return sum(lst)/len(lst)


def make_R_lists(id_list, d, overlay_dict, aggr, intersected_introns):
        s = ""
        aggr_f = {
                "mean": mean,
                "median": median,
        }
        id_list = id_list if not overlay_dict else overlay_dict.keys()
        # Iterate over ids to get bam signal and junctions
        shrinked_introns = dict()
        for k in id_list:
                shrinked_introns_k, shrinked_intronsid = dict(), dict()
                x, y, dons, accs, yd, ya, counts = [], [], [], [], [], [], []
                if not overlay_dict:
                        x, y, dons, accs, yd, ya, counts = d[k]
                        if intersected_introns:
                                x, y = shrink_density(x, y, intersected_introns)
                                shrinked_introns_k, dons, accs = shrink_junctions(dons, accs, intersected_introns)
                                shrinked_introns.update(shrinked_introns_k)
                else:
                        for IDX in overlay_dict[k]:
                                xid, yid, donsid, accsid, ydid, yaid, countsid = d[IDX]
                                if intersected_introns:
                                        xid, yid = shrink_density(xid, yid, intersected_introns)
                                        shrinked_intronsid, donsid, accsid = shrink_junctions(donsid, accsid, intersected_introns)
                                        shrinked_introns.update(shrinked_intronsid)
                                x += xid
                                y += yid
                                dons += donsid
                                accs += accsid
                                yd += ydid
                                ya += yaid
                                counts += countsid
                        if aggr and "_j" not in aggr:
                                x = d[overlay_dict[k][0]][0]
                                y = list(map(aggr_f[aggr], zip(*(d[IDX][1] for IDX in overlay_dict[k]))))
                                if intersected_introns:
                                        x, y = shrink_density(x, y, intersected_introns)
                        #dons, accs, yd, ya, counts = [], [], [], [], []
                s += """
                density_list[["%(IDX)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
                junction_list[["%(IDX)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
                """ %({
                        'IDX': k,
                        'x' : ",".join(map(str, x)),
                        'y' : ",".join(map(str, y)),
                        'dons' : ",".join(map(str, dons)),
                        'accs' : ",".join(map(str, accs)),
                        'yd' : ",".join(map(str, yd)),
                        'ya' : ",".join(map(str, ya)),
                        'counts' : ",".join(map(str, counts))
                })
        
        return s


def plot(R_script):
        p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
        p.communicate(input=R_script.encode('utf-8'))
        p.stdin.close()
        p.wait()
        return


def get_debug_info():
        """
        Return useful debug information:
        - OS info
        - Linux distribution
        - Python version
        - program version
        - output of R `sessionInfo()`
        """
        # get system info
        import platform
        system = platform.system()
        info = OrderedDict()
        info["OS"] = "{}-{}".format(system, platform.machine())
        if system == "Linux":
                release = sp.check_output(["lsb_release", "-ds"])
                info["Distro"] = release.strip().decode('utf-8')
        info["Python"] = platform.python_version()
        # info["ggsashimi"] = __version__
        print(get_version())
        print('')
        maxlen = max(map(len, info.keys()))
        for k,v in info.items():
                print("{:{width}}: {:>}".format(k ,v, width=maxlen))
        print('')

        # get R session info
        r_exec = ';'.join([
                "library(ggplot2)",
                "library(grid)",
                "library(gridExtra)",
                "library(data.table)",
                "library(gtable)",
                "sessionInfo()",
        ])

        r_command = "R --vanilla --slave -e '{}'".format(r_exec)
        r_info = sp.check_output(r_command, shell=True, stderr=sp.STDOUT)
        print(r_info.strip().decode('utf-8'))


def main():

        strand_dict = {"plus": "+", "minus": "-"}

        parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
        
        parser.add_argument("-b", "--bam", type=str, required=True,
                help="""
                Individual bam file or file with a list of bam files.
                In the case of a list of files the format is tsv:
                1col: IDX for bam file,
                2col: path of bam file,
                3+col: additional columns
                """)
        parser.add_argument("-c", "--coordinates", type=str, required=True,
                help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
        parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
                help="Prefix for plot file name [default=%(default)s]")
        parser.add_argument("-S", "--out-strand", type=str, dest="out_strand", default="both",
                help="Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=%(default)s]")
        parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
                help="Minimum number of reads supporting a junction to be drawn [default=1]")
        parser.add_argument( "--min-junction-ratio", type=float, default=0.05, dest="min_fraction",
                help="Minimum number of reads supporting a junction to be drawn [default=0.05]")
        parser.add_argument("-j", "--junctions-bed", type=str, dest = "junctions_bed", default="",
                help="Junction BED file name [default=no junction file]")
        parser.add_argument("-g", "--gtf",
                help="Gtf file with annotation (only exons is enough)")
        parser.add_argument("--mi", dest = "mi",
                help="MI summary file")
        parser.add_argument("-s", "--strand", default="NONE", type=str,
                help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
        parser.add_argument("--shrink", action="store_true",
                help="Shrink the junctions by a factor for nicer display [default=%(default)s]")
        parser.add_argument("-O", "--overlay", type=int,
                help="Index of column with overlay levels (1-based)")
        parser.add_argument("-A", "--aggr", type=str, default="",
                help="""Aggregate function for overlay: <mean> <median> <mean_j> <median_j>.
                        Use mean_j | median_j to keep density overlay but aggregate junction counts [default=no aggregation]""")
        parser.add_argument("-C", "--color-factor", type=int, dest="color_factor",
                help="Index of column with color levels (1-based)", default=3)
        parser.add_argument("--alpha", type=float, default=1.0,
                help="Transparency level for density histogram [default=%(default)s]")
        parser.add_argument("-P", "--palette", type=str,
                help="Color palette file. tsv file with >=1 columns, where the color is the first column. Both R color names and hexadecimal values are valid")
        parser.add_argument("-L", "--labels", type=int, dest="labels", default=3,
                help="Index of column with labels (1-based) [default=%(default)s]")
        parser.add_argument("--fix-y-scale", default=False, action="store_true", dest = "fix_y_scale",
                help="Fix y-scale across individual signal plots [default=%(default)s]")
        parser.add_argument("--height", type=float, default=2,
                help="Height of the individual signal plot in inches [default=%(default)s]")
        parser.add_argument("--ann-height", type=float, default=2.5, dest="ann_height",
                help="Height of annotation plot in inches [default=%(default)s]")
        parser.add_argument("--width", type=float, default=10,
                help="Width of the plot in inches [default=%(default)s]")
        parser.add_argument("--base-size", type=float, default=14, dest="base_size",
                help="Base font size of the plot in pch [default=%(default)s]")
        parser.add_argument("-F", "--out-format", type=str, default="pdf", dest="out_format",
                help="Output file format: <pdf> <svg> <png> <jpeg> <tiff> [default=%(default)s]")
        parser.add_argument("-R", "--out-resolution", type=int, default=300, dest="out_resolution",
                help="Output file resolution in PPI (pixels per inch). Applies only to raster output formats [default=%(default)s]")

        global args

        args = parser.parse_args()


        if args.aggr and not args.overlay:
                print("ERROR: Cannot apply aggregate function if overlay is not selected.")
                exit(1)

        bam_dict, overlay_dict, color_dict, id_list, label_dict = {"+":OrderedDict()}, OrderedDict(), OrderedDict(), [], OrderedDict()

        if args.strand != "NONE": 
                bam_dict["-"] = OrderedDict()
        if args.junctions_bed != "": 
                junctions_list = []

        ALL_transcripts = []

        for IDX, bam, overlay_level, color_level, label_text in read_bam_input(args.bam, args.overlay, args.color_factor, args.labels):
                
                if not os.path.isfile(bam):
                        continue
                
                a, junctions, txs = read_bam(bam, args.coordinates, args.strand) 
                ALL_transcripts.extend(txs)

                if a.keys() == ["+"] and all(map(lambda x: x==0, list(a.values()[0]))):
                        print("WARN: Sample {} has no reads in the specified area.".format(IDX))
                        continue
                id_list.append(IDX)
                label_dict[IDX] = label_text
                for strand in a:
                        # Store junction information
                        if args.strand == "NONE" or args.out_strand == 'both' or strand == strand_dict[args.out_strand]:
                                if args.junctions_bed:
                                        max_v = 0
                                        for k, v in zip(junctions[strand].keys(), junctions[strand].values()):
                                                if v > max_v: max_v = v
                                        
                                        for k, v in zip(junctions[strand].keys(), junctions[strand].values()):
                                                if v >= args.min_coverage and float(v/max_v) >= args.min_fraction:
                                                        _chrom = args.coordinates.split(':')[0]
                                                        jxn_string = f'{_chrom}\t{k[0]}\t{k[1]}\t{IDX}\t{v}\t{strand}'
                                                        junctions_list.append(jxn_string)

                        bam_dict[strand][IDX] = prepare_for_R(a[strand], junctions[strand], args.coordinates, args.min_coverage, args.min_fraction)
                if color_level is None:
                        color_dict.setdefault(IDX, IDX)
                if overlay_level is not None:
                        overlay_dict.setdefault(overlay_level, []).append(IDX)
                        label_dict[overlay_level] = overlay_level
                        color_dict.setdefault(overlay_level, overlay_level)
                if overlay_level is None:
                        color_dict.setdefault(IDX, color_level)

        # No bam files
        if not bam_dict["+"]:
                print("ERROR: No available bam files.")
                exit(1)

        # Write junctions to BED
        if args.junctions_bed:
                if not args.junctions_bed.endswith('.bed'):
                        args.junctions_bed = args.junctions_bed + '.bed'
                jbed = open(args.junctions_bed, 'w')
                jbed.write('\n'.join(sorted(junctions_list)))
                jbed.close()

        if args.gtf:
                transcripts, exons = read_gtf(args.gtf, args.coordinates, ALL_transcripts)

        if args.out_format not in ('pdf', 'png', 'svg', 'tiff', 'jpeg'):
                print("ERROR: Provided output format '%s' is not available. Please select among 'pdf', 'png', 'svg', 'tiff' or 'jpeg'" % args.out_format)
                exit(1)

        # Iterate for plus and minus strand
        for strand in bam_dict:

                # Output file name (allow tiff/tif and jpeg/jpg extensions)
                if args.out_prefix.endswith(('.pdf', '.png', '.svg', '.tiff', '.tif', '.jpeg', '.jpg')):
                        out_split = os.path.splitext(args.out_prefix)
                        if (args.out_format == out_split[1][1:] or
                        args.out_format == 'tiff' and out_split[1] in ('.tiff','.tif') or
                        args.out_format == 'jpeg' and out_split[1] in ('.jpeg','.jpg')):
                                args.out_prefix = out_split[0]
                                out_suffix = out_split[1][1:]
                        else:
                                out_suffix = args.out_format
                else:
                        out_suffix = args.out_format
                out_prefix = args.out_prefix + "_" + strand
                if args.strand == "NONE":
                        out_prefix = args.out_prefix
                else:
                        if args.out_strand != "both" and strand != strand_dict[args.out_strand]:
                                continue

                # Find set of junctions to perform shrink
                intersected_introns = None
                if args.shrink:
                        introns = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
                        intersected_introns = list(intersect_introns(introns))


                # *** PLOT *** Define plot height

                R_script = setup_R_script(label_dict)

                # *** PLOT *** Prepare annotation plot only for the first bam file
                arrow_bins = 50
                if args.gtf:
                        # Make introns from annotation (they are shrunk if required)
                        annotation = make_introns(transcripts, exons, intersected_introns)
                        x = list(bam_dict[strand].values())[0][0]
                        if args.shrink:
                                x, _ = shrink_density(x, x, intersected_introns)
                        R_script += gtf_for_ggplot(annotation, args.coordinates,  args.mi, x[0], x[-1], arrow_bins)

                R_script += make_R_lists(id_list, bam_dict[strand], overlay_dict, args.aggr, intersected_introns)

                R_script += """

                plot_list <- list()
                palette = brewer.pal(2, "Dark2") # track colors

                for (bam_index in 1:length(density_list)) {

                        IDX = names(density_list)[bam_index]
                        d = data.table(density_list[[IDX]])
                        
                        dat = d[d$y > 0,]
                        junctions = data.table(junction_list[[IDX]])

                        gp = ggplot(dat) 
                        gp = gp + geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=palette[bam_index], color=palette[bam_index])
                        gp = gp + labs(y=labels[[IDX]], x = '')
                        gp = gp + scale_x_continuous(expand=c(0, 0.25), limits = c(min(d$x), max(d$x)))                        

                        # Aggregate junction counts
                        row_i = c()
                        if (nrow(junctions) >0 ) {

                                junctions$jlabel = as.character(junctions$count)
                                junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),paste(jlabel,collapse=",")), keyby=.(x,xend)], names(junctions))
                                if ("%(args.aggr)s" != "") {
                                        junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),round(%(args.aggr)s(count))), keyby=.(x,xend)], names(junctions))
                                }
                                # The number of rows (unique junctions per bam) has to be calculated after aggregation
                                row_i = 1:nrow(junctions)
                        }


                        for (i in row_i) {

                                j_tot_counts = sum(junctions[['count']])

                                j = as.numeric(junctions[i,1:5])

                                if ("%(args.aggr)s" != "") {
                                        j[3] = ifelse(length(d[x==j[1]-1,y])==0, 0, max(as.numeric(d[x==j[1]-1,y])))
                                        j[4] = ifelse(length(d[x==j[2]+1,y])==0, 0, max(as.numeric(d[x==j[2]+1,y])))
                                }

                                # Find intron midpoint
                                xmid = round(mean(j[1:2]), 1)
                                set.seed(mean(j[3:4]))
                                ymid = max(j[3:4]) * runif(1, 1.2, 1.5)

                                # Thickness of the arch
                                thickness.arch <- max(0.5, (j[5]*3)/j_tot_counts)

                                curve_par = gpar(lwd=0.5, col=palette[bam_index], alpha = 0.6)

                                # Arc grobs

                                nss = i
                                if (nss%%%%2 == 0) {  #bottom
                                        ymid = -runif(1, 0.2, 0.4) * max(j[3:4])
                                        # Draw the arcs
                                        # Left
                                        curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, j[1], xmid, 0, ymid)
                                        # Right
                                        curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, xmid, j[2], 0, ymid)
                                }

                                if (nss%%%%2 != 0) {  #top
                                        # Draw the arcs
                                        # Left
                                        curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, j[1], xmid, j[3], ymid)
                                        # Right
                                        curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                                        gp = gp + annotation_custom(grob = curve, xmid, j[2], j[4], ymid)
                                }

                                # Add junction labels
                                y_off = max(d$y) * 0.2
                                gp = gp + annotate("label", x = xmid, y = (ymid + (ymid/abs(ymid))*y_off), label = as.character(junctions[i,6]),
                                        vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"),
                                        label.size=NA, size = 2.5)

                                gp = gp + theme_bw() + theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

                        }
                        plot_list[[bam_index]] = gp
                }

                plot_list[[(length(density_list) + 1)]] = gtfp                
                argrobs <- wrap_plots(plot_list, ncol = 1)
                ggsave("%(out)s", plot = argrobs, device = "%(out_format)s", width = 9, height = 6, units = "in", dpi = %(out_resolution)s, limitsize = FALSE)
                dev.log = dev.off()

                """ %({
                        "out": "%s.%s" % (out_prefix, out_suffix),
                        "out_format": args.out_format,
                        "out_resolution": args.out_resolution,
                        "args.gtf": float(bool(args.gtf)),
                        "args.aggr": args.aggr.rstrip("_j"),
                        "signal_height": args.height,
                        "ann_height": args.ann_height,
                        "alpha": args.alpha,
                        "fix_y_scale": ("TRUE" if args.fix_y_scale else "FALSE")
                        })

                #if os.getenv('GGSASHIMI_DEBUG') is not None:
                with open("R_script", 'w') as r:
                        r.write(R_script)
                #else:
                plot(R_script)
        exit()




if __name__ == "__main__":
        main()
