# **isoLASER**
[![](https://img.shields.io/badge/isoLASER-v0.0.0.1-blue)](https://test.pypi.org/project/isoLASER/)

[github](https://github.com/gxiaolab/isoLASER/)
_______________________________________

## **About**


isoLASER performs gene-level variant calls, phasing and splicing linkage analysis using third generation RNA sequencing data.

## **Table of contents**
- [Outline](#Outline)
- [Installation](#Installation)
- [Preprocessing](#Preprocessing)
  - [Annotate bam file](#Annotate bam file)
  - [Extract exonic parts from GTF](#Extract exonic parts from GTF)
- [Run isoLASER](#Run isoLASER)
- [Run isoLASER joint](#Run isoLASER joint)
- [Make a nigiri plot!](#Make a nigiri plot!)
- [Output](#Output)
- [Debug](#Debug)



## **Outline**


## **Installation**

isoLASER is available through **PyPi**. To download simply type:

```
pip install isoLASER
```

The download was tested with PyPi version >= 20.0.1.

Alternatively, you can clone this **GitHub** repository:

```
git clone git@github.com:gxiaolab/isoLASER.git 
cd isoLASER
python -m build
pip install .
```

You can also download the **Singularity** container:

```
singularity pull library://giovas/collection/s6
singularity exec s5_latest.sif isoLASER
```

If successful, the program is ready to use. The installation incorporates console script entry points to directly call isoLASER:

```
isoLASER --help
```
_______________________________________



## *Preprocessing* 

Long-read RNA sequencing is notorious for its high base-calling error rate. As such, it is important to clean and preprocess the data to discard false transcripts resulting from misalignment, bad consensus, truncation, and other technical artifacts.   

isoLASER a GTF file as input, ideally built using a long read annotation software such as Talon, Clair, Bambu, Espresso, or similar. 

### Talon pipeline

Talon + Transcript clean 
For more details see (link)[www.=]

```
python TranscriptClean-2.0.4/accessory_scripts/get_SJs_from_gtf.py \
                        --f {input.gtf} \
                        --g {input.ref} \
                        --o {output.jxn}

python TranscriptClean-2.0.4/TranscriptClean.py \
                        --correctMismatches False \
                        --correctIndels False \
                        --threads {params.threads} \
                        --sam {output.sam} \
                        --genome {input.ref} \
                        --spliceJns {input.jxn} \
                        --outprefix {params.prefix} \
                        --tmpDir {params.tmpdir}

talon_label_reads \
                        --f {input.sam} \
                        --g {input.ref} \
                        --o {params.prefix} \
                        --tmpDir {params.tmp}

talon_initialize_database \
                        --f talon/$chrom.gtf \
                        --g {params.genome_assembly} \
                        --a {params.annot_version} \
                        --idprefix {params.dataset_id} \
                        --o talon/$chrom

talon \
                        --f {output.Config} \
                        --db {input.db} \
                        --build {params.genome_assembly} \
                        --threads 24 \
                        --tmpDir {params.tmpdir} \
                        --o {params.prefix}

 talon_filter_transcripts \
                        --db {input.db} \
                        -a {params.dataset_id} \
                        --includeAnnot \
                        --minCount 5 \
                        --minDatasets 1 \
                        --o {output.whitelist}

  talon_create_GTF 
                        --db {input.db} \
                        -a {params.annot_version} \
                        --build {params.genome_assembly} \
                        --whitelist {input.whitelist} \
                        --observed \
                        --o {params.prefix}

```


### Annotate bam file

Obtain a reference transcriptome file to align the reads to. 
This step serves to assign transcript ids to every read of the target bam file

```
isolaser_convert_gtf_to_fasta
```

Align agains the newly generated reference

```
minimap2
```

Filter for secondary, supplementary and trans-gene reads whilst annotating with transcript ids. 
Transcript ids and gene names are saved in the `ZG` and `ZT` tags.

```
isolaser_filter_and_annotate
```

### Extract exonic parts from GTF

isoLASER uses a exon-centric approach to analyze splicing and exonic-parts are great granular appraoch to understand local splicing changes. 

```
isolaser_extract_exonic_parts
```


## *Run isoLASER*

```
isoLASER -b <file.bam> -o <output_prefix> -t <transcriptome.db> -f <reference.fa>
```

The output is very extensive and includes information that is only relevant for the joint analysis or plotting. 
To obtain the significant allele-speicifc events (cis-drected splicing events) use the filter function: 

```
isolaser_filter -m <output_prefix.mi_summary.tsv> -o <output_prefix.mi_summary.filtered.tsv>
```

## *Run isoLASER joint*

Wrapper of GATK functions to merge the variant calls from different samples 

```
isoLASER_combine_vcf -f <fofn.tsv> -o <output_prefix> 
```
Perform joint analysis
```
isoLASER_joint -f <fofn.tsv> -o <output_prefix> -t <transcriptome.db> 
```

## *Make a nigiri plot!*

```
nigiri
```

## **Output**


## **Debug**


