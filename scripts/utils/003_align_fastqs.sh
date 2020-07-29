

# checking if the required directory structure exists
if [[ ! -e data ]]; then mkdir -p data; fi
if [[ ! -e data/fastq ]]; then mkdir -p data/fastq; fi