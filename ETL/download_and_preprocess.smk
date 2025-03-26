# Note: To avoid the ILP CBC error, run Snakemake with the greedy scheduler:
#       snakemake --scheduler greedy -j <jobs>
SAMPLES = ["SRR5088815", "SRR5088816"]

rule all:
    input:
        expand("data/{sample}/{sample}_cdr3.out", sample=SAMPLES)

rule prefetch:
    output:
        sra="data/{sample}/{sample}.sra"
    threads: 1
    shell:
        """
        docker run --rm -v /Users/talshor/Projects/Numenos_ex/data:/data srr_download_data \
            prefetch {wildcards.sample} --output-directory /data
        """

rule fastq_dump:
    input:
        sra="data/{sample}/{sample}.sra"
    output:
        read1="data/{sample}/{sample}_1.fastq",
        read2="data/{sample}/{sample}_2.fastq"
    threads: 1
    shell:
        """
        docker run --rm -v /Users/talshor/Projects/Numenos_ex/data:/data srr_download_data \
            fastq-dump --split-files {wildcards.sample} --outdir /data/{wildcards.sample}
        """

rule trust4:
    input:
        read1="data/{sample}/{sample}_1.fastq",
        read2="data/{sample}/{sample}_2.fastq"
    output:
        cdr3="data/{sample}/{sample}_cdr3.out",
        report="data/{sample}/{sample}_report.tsv",
        toassemble1="data/{sample}/{sample}_toassemble_1.fq",
        toassemble2="data/{sample}/{sample}_toassemble_2.fq",
        airr="data/{sample}/{sample}_airr.tsv",
        airr_align="data/{sample}/{sample}_airr_align.tsv",
        annot="data/{sample}/{sample}_annot.fa",
        assembled="data/{sample}/{sample}_assembled_reads.fa",
        raw="data/{sample}/{sample}_raw.out",
        final="data/{sample}/{sample}_final.out"
    threads: 4
    shell:
        """
        docker run --rm -v /Users/talshor/Projects/Numenos_ex/data/{wildcards.sample}:/data trust4_docker \
            run-trust4 -f /reference/hg38_bcrtcr.fa --ref /reference/human_IMGT+C.fa \
            -1 data/{wildcards.sample}_1.fastq -2 data/{wildcards.sample}_2.fastq -o /data/{wildcards.sample}
        """
