readme: README.rst

README.rst: README.md
	pandoc -f markdown -t rst README.md > README.rst

test: pooled unpooled dirty missing

pooled: tests/pooled/134501_LibID134601_GAP_BRF_H5TT7DRX3_GCTTGTCA-GTATGTTC_S82_L002_r1.fastq.gz
tests/pooled/134501_LibID134601_GAP_BRF_H5TT7DRX3_GCTTGTCA-GTATGTTC_S82_L002_r1.fastq.gz:
	mkdir -p tests/pooled && \
	tcdemux \
	--sample_data data/samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/pooled \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \
	2>&1 | tee tests/pooled/sm.log

unpooled: tests/unpooled/134486_LibID134586_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S10_L002_r1.fastq.gz
tests/unpooled/134486_LibID134586_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S10_L002_r1.fastq.gz:
	mkdir -p tests/unpooled && \
	tcdemux \
	--sample_data data/unpooled_samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/unpooled \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16
	2>&1 | tee tests/unpooled/sm.log

dirty: tests/dirty/sm.log
tests/dirty/sm.log:
	mkdir -p tests/dirty && \
	tcdemux \
	--sample_data data/dirty_samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/dirty \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \
	2>&1 | tee tests/dirty/sm.log

missing: tests/missing/this_sample_exists.unpaired.fastq.gz
tests/missing/this_sample_exists.unpaired.fastq.gz:
	mkdir -p tests/missing && \
	tcdemux \
	--sample_data data/samples_missing.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/missing \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \
	2>&1 | tee tests/missing/sm.log

graph: assets/graph.svg assets/external_only_graph.svg

assets/graph.svg:
	snakemake -n \
	-s tcdemux/Snakefile \
	--cores 12 \
	--rulegraph \
	--forceall \
	--config \
	sample_data_file=data/samples_graph.csv \
	read_directory=data/test_data \
	adaptor_files=data/adaptors/TruSeq3-PE-2.fa \
	outdir=tests/graph \
	mem_gb=12 \
	threads=12 \
	keep_intermediate_files=False \
	| grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
	| sed 's/label = "7"/label = "parse_cutadapt_stats"/g' \
	| sed 's/label = "8"/label = "cutadapt"/g' \
	| sed 's/label = "9"/label = "write_barcode_file"/g' \
	| dot -Tsvg \
	> assets/graph.svg

assets/external_only_graph.svg:
	snakemake -n \
	-s tcdemux/Snakefile \
	--cores 12 \
	--rulegraph \
	--forceall \
	--config \
	sample_data_file=data/unpooled_samples_graph.csv \
	read_directory=data/test_data \
	adaptor_files=data/adaptors/TruSeq3-PE-2.fa \
	outdir=tests/graph \
	mem_gb=12 \
	threads=12 \
	keep_intermediate_files=False \
	| grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
	| dot -Tsvg \
	> assets/external_only_graph.svg
