docker build -t srr_download_data .

docker run --rm -v $(pwd)/data:/data srr_download_data \
    SRR5088813 --output-directory /data