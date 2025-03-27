SAMPLES = ["SRR5088829", "SRR5088830", "SRR5088831", "SRR5088832"]

rule all:
    input:
        expand("data/{sample}/{sample}_cdr3.out", sample=SAMPLES),
        #Wexpand("data/{sample}/{sample}_fastqc_raw.html", sample=SAMPLES),
        
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

rule fastqc_raw:
    input:
        read1="data/{sample}/{sample}_1.fastq",
        read2="data/{sample}/{sample}_2.fastq"
    output:
        report="data/{sample}/{sample}_fastqc_raw.html",
        zip="data/{sample}/{sample}_fastqc_raw.zip"
    threads: 2
    shell:
        """
        docker run --rm -v /Users/talshor/Projects/Numenos_ex/data:/data fastqc_docker \
            fastqc -o /data/{wildcards.sample} --threads {threads} \
            data/{wildcards.sample}/{wildcards.sample}_1.fastq \
            data/{wildcards.sample}/{wildcards.sample}_2.fastq
        """

# The quality is really good - no need for this
# rule fastp:
#     input:
#         read1="data/{sample}/{sample}_1.fastq",
#         read2="data/{sample}/{sample}_2.fastq"
#     output:
#         read1="data/{sample}/{sample}_1.trimmed.fastq",
#         read2="data/{sample}/{sample}_2.trimmed.fastq",
#     threads: 2
#     shell:
#         """
#         docker run --rm -v /Users/talshor/Projects/Numenos_ex/data:/data fastp_docker \
#             fastp -i data/{wildcards.sample}/{wildcards.sample}_1.fastq \
#             -I data/{wildcards.sample}/{wildcards.sample}_2.fastq \
#             -o data/{wildcards.sample}/{wildcards.sample}_1.trimmed.fastq \
#             -O data/{wildcards.sample}/{wildcards.sample}_2.trimmed.fastq \
#             --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
#             --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --length_required 50
#         """

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
    threads: 2
    shell:
        """
        docker run --rm -v /Users/talshor/Projects/Numenos_ex/data/{wildcards.sample}:/data trust4_docker \
            run-trust4 -f /reference/hg38_bcrtcr.fa --ref /reference/human_IMGT+C.fa \
            -1 data/{wildcards.sample}_1.fastq -2 data/{wildcards.sample}_2.fastq -o /data/{wildcards.sample} \
            -t {threads}
        """
