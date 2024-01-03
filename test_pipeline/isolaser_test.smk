# Snakefile

configfile: "config.yaml"


rule all:
    input:
        expand("out/isolaser.dlpfc_{i}.mi_summary.filtered.tab", i=[1, 2, 3]) + 
        ["out/isolaser_joint.dlpfc.merged.mi_summary.tab"]


rule convert_gff_to_fasta:
    input:
        gff = config["GTF"],
        ref = config["REF"],
    output:
        transcriptome = "data/talon_observedOnly.sorted.fa"
    shell:
        """
        mkdir -p data

        isolaser_convert_gtf_to_fasta \
            -g {input.gff} \
            -f {input.ref} \
            -o {output.transcriptome}
        """


rule extract_exonic_parts:
    params:
        db = "data/talon_observedOnly_db"
    input:
        gff = config["GTF"],
    output:
        transcriptome = "data/talon_observedOnly_db/GeneList.tsv"
    shell:
        """
        isolaser_extract_exon_parts \
            -g {input.gff} \
            -o {params.db}
        """


rule map_to_transcriptome:
    input:
        bam = "bam/dlpfc_{i}.bam",
        transcriptome = "data/talon_observedOnly.sorted.fa"
    output:
        fq  = "bam/dlpfc_{i}.fq",
        sam = "bam/dlpfc_{i}.transcriptome.sam"
    shell:
        """
        samtools fastq {input.bam} > {output.fq}

        minimap2 -t 16 -ax splice:hq \
            -uf --MD {input.transcriptome} {output.fq} > {output.sam}
        """


rule annotate_reads:
    input:
        bam = "bam/dlpfc_{i}.bam",
        sam = "bam/dlpfc_{i}.transcriptome.sam",
        gff = config["GTF"]
    output:
        bam = "bam/dlpfc_{i}.annotated.bam"
    shell:
        """
        isolaser_annotate \
            -b {input.bam} \
            -t {input.sam} \
            -g {input.gff} \
            -o {output.bam}
        """


rule bam_sort_index:
    input:
        bam = "bam/dlpfc_{i}.annotated.bam"
    output:
        bam = "bam/dlpfc_{i}.annotated.sorted.bam"
    shell:
        """
        samtools sort {input.bam} -o {output.bam}

        samtools index {output.bam}
        """


rule isolaser_run:
    params:
        prefix = "out/isolaser.dlpfc_{i}",
        transcriptome_db = "data/talon_observedOnly_db",
    input:
        bam = "bam/dlpfc_{i}.annotated.sorted.bam",       
        transcriptome_db_lst = "data/talon_observedOnly_db/GeneList.tsv",
        ref = config["REF"]
    output:
        mi_summary = "out/isolaser.dlpfc_{i}.mi_summary.tab.gz"
    shell:
        """
        isolaser \
            -b {input.bam} \
            -o {params.prefix} \
            -t {params.transcriptome_db} \
            -f {input.ref} &> log/isolaser.dlpfc_{wildcards.i}.log

        less {params.prefix}.mi_summary.tab | sort -nk2 -k1,1V -k2,2n -k3,3n | bgzip -c -f > {output.mi_summary}

        tabix -f -p bed {output.mi_summary}
        """


rule isolaser_filter:
    input:
        mi = "out/isolaser.dlpfc_{i}.mi_summary.tab.gz"
    output:
        mi = "out/isolaser.dlpfc_{i}.mi_summary.filtered.tab"
    shell:
        """
        isolaser_filter \
            -m {input.mi} \
            -o {output.mi}
        """


rule write_fofn:
    input:
        mi = expand("out/isolaser.dlpfc_{i}.mi_summary.tab.gz", i=[1, 2, 3])
    output:
        fofn = "out/fofn.txt"
    shell:
        """
        echo -n > {output.fofn}
        for i in 1 2 3;
        do
            bamf="bam/dlpfc_${{i}}.annotated.sorted.bam"
            mif="out/isolaser.dlpfc_${{i}}.mi_summary.tab.gz"
            vcf="out/isolaser.dlpfc_${{i}}.gvcf"

            printf "${{bamf}}\t${{mif}}\t${{vcf}}\n" >> {output.fofn}
        done
        cat {output.fofn}    
        """


rule isolaser_combine_vcf:
    params:
        prefix = "out/isolaser_joint.dlpfc"
    input:
        fofn = "out/fofn.txt",
        ref = config["REF"],
    output:
        vcf = "out/isolaser_joint.dlpfc.merged.genotyped.gvcf.gz",
    shell:
        """ 
        isolaser_combine_vcf \
            -i {input.fofn} \
            -f {input.ref} \
            -o {params.prefix} \
        """


rule isolaser_joint:
    params:
        transcriptome_db = "data/talon_observedOnly_db",
    input:
        fofn = "out/fofn.txt",
        vcf  = "out/isolaser_joint.dlpfc.merged.genotyped.gvcf.gz",
        transcriptome_db_lst = "data/talon_observedOnly_db/GeneList.tsv",
    output:
        mi = "out/isolaser_joint.dlpfc.merged.mi_summary.tab"
    shell:
        """ 
        isolaser_joint \
            -i {input.fofn} \
            -t {params.transcriptome_db} \
            -v {input.vcf} \
            -o {output.mi} 
        """
