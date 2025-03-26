# docker build -t srr_download_data .

docker run --rm -v $(pwd)/../data:/data srr_download_data \
    prefetch SRR5088813 --output-directory /data

docker run --rm -v $(pwd)/../data:/data srr_download_data \
    fastq-dump --split-files SRR5088813 --outdir /data/SRR5088813

docker run --rm -v $(pwd)/../data/SRR5088813:/data trust4_docker \
    run-trust4 -f /reference/hg38_bcrtcr.fa --ref /reference/human_IMGT+C.fa -1 data/SRR5088813_1.fastq -2 data/SRR5088813_2.fastq -o data/SRR5088813