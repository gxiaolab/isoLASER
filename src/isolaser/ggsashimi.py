#!/usr/bin/env python
import subprocess as sp
import sys, re, copy, os, codecs, gzip
from argparse import ArgumentParser, Action as ArgParseAction
from collections import OrderedDict, Counter, defaultdict
import pysam
import numpy as np
import HTSeq


def parse_coordinates(c):
        c = c.replace(",", "")
        chrom = c.split(":")[0]
        start, end = c.split(":")[1].split("-")
        start, end = int(start) - 1, int(end)
        return chrom, start, end



def count_operator(CIGAR_op, CIGAR_len, pos, start, end, cov, junctions):
        
        # Match
        if CIGAR_op in [0, 7, 8]: 
                for i in range(pos, pos + CIGAR_len):
                        if i < start or i >= end:
                                continue
                        ind = i - start
                        cov[ind] += 1

        # Insertion or Soft-clip
        if CIGAR_op == 1 or CIGAR_op == 4:
                return pos

        # Deletion
        if CIGAR_op == 2:
                pass

        # Junction
        if CIGAR_op == 3:
                don = pos
                acc = pos + CIGAR_len
                if don > start and acc < end:
                        junctions[(don,acc)] = junctions.setdefault((don,acc), 0) + 1

        pos = pos + CIGAR_len

        return pos


def read_bam(bam_file, coords):

        chrom, start, end = parse_coordinates(coords)

        # Initialize coverage array and junction dict
        cov = {"+" : [0] * (end - start)}
        junctions = {"+": OrderedDict()}
        

        samfile = pysam.AlignmentFile(bam_file)
        transcripts = []

        for read in samfile.fetch(chrom, start, end):
                
                # Move forward if read is unmapped
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue           

                if any(map(lambda x: x in read.cigarstring, ["H", "P"])):
                        continue

                read_strand = "+"
               
                pos = read.reference_start + 1 
                for match_str, match_len in read.cigar:
                        pos = count_operator(match_str, match_len, pos, start, end, cov[read_strand], junctions[read_strand])
               
                transcripts.append(read.get_tag("ZT"))

        samfile.close()
        return cov, junctions, transcripts



def read_bam_input(f):
        with codecs.open(f, encoding='utf-8') as openf:
                for line in openf:
                        line_sp = line.strip().split("\t")
                        assert len(line_sp) == 3, "ERROR: The input file should have 3 columns: id, bam, label"
                        sm_id, bam, label = line_sp
                        assert os.path.isfile(bam), f"ERROR: The bam file does not exist {bam}."
                        yield sm_id, bam, label


def prepare_for_R(a, junctions, c, min_cov, min_freq):

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

        for (don, acc), jxn_count in junctions.items():

                # Do not add junctions with less than defined coverage
                if jxn_count < min_cov or (jxn_count/max_n < min_freq):
                        continue

                dons.append(don)
                accs.append(acc)
                counts.append(jxn_count)

                yd.append( a[ don - start -1 ])
                ya.append( a[ acc - start +1 ])

        return x, y, dons, accs, yd, ya, counts


def read_gtf(GTF_file, c, tx_id_list):

        n = len(tx_id_list)
        transcript_frequency = {}

        for k, v in Counter(tx_id_list).items():
                transcript_frequency[f'"{k}"'] = float(v/n)

        exons = OrderedDict()
        transcripts = OrderedDict()

        chrom, start, end = parse_coordinates(c)
        c = HTSeq.GenomicInterval(chrom, start, end - 1, ".")

        tabixfile = pysam.TabixFile(GTF_file, parser = pysam.asGTF())
        for record in tabixfile.fetch(chrom, start, end):
                print(record.feature, record.transcript_id )
                if record.feature not in ("transcript", "exon"):
                        continue
                try:
                        transcript_id = record.transcript_id
                except KeyError:
                        print("ERROR: 'transcript_id' attribute is missing in the GTF file.")
                        exit(1)

                if transcript_id not in transcript_frequency.keys():
                        continue

                print(record) ; break

                record_iv = HTSeq.GenomicInterval(record.contig, record.start, record.end, record.strand)
                
                if transcript_frequency[transcript_id] >= 0.05 or "ENST" in transcript_id:

                        strand = '"' + strand + '"'
                        if record.feature == "transcript":
                                if record_iv.overlaps(c):
                                        transcripts[transcript_id] = record_iv
                        if record.feature == "exon":
                                if c.contains(record_iv):
                                        exons.setdefault(transcript_id, []).append((record_iv))

        return transcripts, exons


def make_introns(transcripts, exons, intersected_introns=None):
        new_transcripts = copy.deepcopy(transcripts)
        new_exons = copy.deepcopy(exons)
        introns = OrderedDict()
       
        for tx, tx_iv in new_transcripts.items():
                intron_start = tx_iv.start
                ex_end = 0
                for ex_iv in sorted(new_exons.get(tx, [])):
                        intron_end = ex_iv.start
                        if tx_iv.start < ex_iv.start:
                                intron_vi = HTSeq.GenomicInterval(tx_iv.chrom, intron_start, intron_end, tx_iv.strand)
                                introns.setdefault(tx, []).append(intron_vi)
                        intron_start = ex_end
                if tx_iv.end > ex_iv.end:
                        intron_vi = HTSeq.GenomicInterval(tx_iv.chrom, intron_start, tx_iv.end, tx_iv.strand)
                        introns.setdefault(tx, []).append(intron_vi)
        d = {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}
        return d


def gtf_for_ggplot(annotation, coords, mi_file, arrow_bins):
        chrom, start, end = parse_coordinates(coords)
        arrow_space = int((end - start)/arrow_bins)
        

        s = """
        # data table with exons
        ann_list = list(
                "exons" = data.table(),
                "introns" = data.table()
        )
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
                "exon_s" : ",".join(map(str, (v.start for vs in annotation["exons"].values() for v in vs))),
                "exon_e" : ",".join(map(str, (v.end for vs in annotation["exons"].values() for v in vs))),
                "strand" : ",".join(map(str, (v.strand for vs in annotation["exons"].values() for v in vs))),
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
                        "intron_s" : ",".join(map(str, (v.start for vs in annotation["introns"].values() for v in vs))),
                        "intron_e" : ",".join(map(str, (v.end for vs in annotation["introns"].values() for v in vs))),
                        "strand"   : ",".join(map(str, (v.strand for vs in annotation["introns"].values() for v in vs))),
                        "arrow_space" : 5000,
                })

        s += """
        gtfp = ggplot()

        head(ann_list[['exons']])
        head(ann_list[['introns']])
        
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
        gtfp = gtfp + labs(y=NULL, x=%s)
        gtfp = gtfp + theme_bw() 
        """ %(start, end, chrom)

        return s


def setup_R_script():
        s = """
        library(ggplot2)
        library(grid)
        library(gridExtra)
        library(data.table)
        library(gtable)
        library(patchwork)
        library(RColorBrewer)

        density_list = list()
        junction_list = list()
        """ 
        return s



def make_R_lists(bam_dict):
        s = ""
        for sample_id, lists in bam_dict.items():              
                xid, yid, donsid, accsid, ydid, yaid, countsid = lists
                
                s += """
                density_list[["%(IDX)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
                junction_list[["%(IDX)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))

                """ %({
                        'IDX': sample_id,
                        'x' : ",".join(map(str, xid)),
                        'y' : ",".join(map(str, yid)),
                        'dons' : ",".join(map(str, donsid)),
                        'accs' : ",".join(map(str, accsid)),
                        'yd' : ",".join(map(str, ydid)),
                        'ya' : ",".join(map(str, yaid)),
                        'counts' : ",".join(map(str, countsid))
                })
        
        return s


def plot(R_script):
        p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
        p.communicate(input=R_script.encode('utf-8'))
        p.stdin.close()
        p.wait()
        return




def main():
        parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
        
        parser.add_argument("-b", "--bam", type=str, required=True,
                help="tsv file with columns: id,bam_file,label")
        parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
                help="Prefix for plot file name [default=%(default)s]")
        parser.add_argument("-c", "--coordinates", type=str, required=True,
                help="Genomic coordinates in the format chr:start-end")
        parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
                help="Minimum number of reads supporting a junction to be drawn [default=1]")
        parser.add_argument( "--min-junction-ratio", type=float, default=0.05, dest="min_fraction",
                help="Minimum number of reads supporting a junction to be drawn [default=0.05]")
        parser.add_argument("-g", "--gtf",
                help="Gtf file with annotation (only exons is enough)")
        parser.add_argument("--mi", dest = "mi",
                help="MI summary file")
        
        
        global args

        args = parser.parse_args()

        bam_dict = {"+" : OrderedDict()}

        ALL_transcripts = []

        for IDX, bam, label_text in read_bam_input(args.bam):                
                a, junctions, txs = read_bam(bam, args.coordinates) 
                ALL_transcripts.extend(txs)
                for strand in a:
                        bam_dict[strand][label_text] = prepare_for_R(a[strand], junctions[strand], args.coordinates, args.min_coverage, args.min_fraction)
     

        if args.gtf:
                transcripts, exons = read_gtf(args.gtf, args.coordinates, ALL_transcripts)
        
        # Iterate for plus and minus strand
        for strand in bam_dict:

                R_script = setup_R_script()

                # *** PLOT *** Prepare annotation plot only for the first bam file
                if args.gtf:
                        annotation = make_introns(transcripts, exons)
                        R_script += gtf_for_ggplot(annotation, args.coordinates, args.mi, arrow_bins = 50)

                R_script += make_R_lists(bam_dict[strand])

                R_script += """

                plot_list <- list()
                palette = brewer.pal(3, "Dark2") # track colors

                for (bam_index in 1:length(density_list)) {

                        IDX = names(density_list)[bam_index]
                        d = data.table(density_list[[IDX]])
                        dat = d[d$y > 0,]

                        gp = ggplot(dat) 
                        gp = gp + geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=palette[bam_index], color=palette[bam_index])
                        gp = gp + labs(y=IDX, x = NULL)
                        gp = gp + scale_x_continuous(expand=c(0, 0.25), limits = c(min(d$x), max(d$x)))                        

                        # Aggregate junction counts
                        junctions = data.table(junction_list[[IDX]])

                        row_i = c()
                        if (nrow(junctions) >0 ) {
                                junctions$jlabel = as.character(junctions$count)
                                junctions = setNames(junctions[,.(max(y), max(yend),round(median(count)),paste(jlabel,collapse=",")), keyby=.(x,xend)], names(junctions))
                                # The number of rows (unique junctions per bam) has to be calculated after aggregation
                                row_i = 1:nrow(junctions)
                        }

                        for (i in row_i) {

                                j_tot_counts = sum(junctions[['count']]) 
                                j = as.numeric(junctions[i,1:5])

                                # Find intron midpoint
                                xmid = round(mean(j[1:2]), 1)
                                set.seed(1)
                                ymid = max(j[3:4]) * runif(1, 1.2, 1.5)

                                curve_par = gpar(lwd=0.75, col=palette[bam_index], alpha = 0.6)

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
                                gp = gp + annotate("label", x = xmid, y = (ymid + (ymid/abs(ymid))*y_off), 
                                        label = as.character(junctions[i,6]),
                                        vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"),
                                        label.size=NA, size = 2.5)

                                gp = gp + theme_bw() + theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

                        }
                        plot_list[[bam_index]] = gp
                }

                plot_list[[(length(density_list) + 1)]] = gtfp                
                
                pdf('out/sea.pdf', width = 9, height = 6)
                print(plot_list[[1]])
                print(plot_list[[2]])
                print(plot_list[[3]])
                dev.off()

                argrobs <- wrap_plots(plot_list, ncol = 1)
                ggsave("%(out)s", plot = argrobs, device = "pdf", width = 9, height = 6, units = "in")
                dev.log = dev.off()

                """ %({"out": "%s" % (args.out_prefix)})

                #if os.getenv('GGSASHIMI_DEBUG') is not None:
                with open("R_script", 'w') as r:
                        r.write(R_script)
                #else:
                plot(R_script)
        exit()




if __name__ == "__main__":
        main()
