SHELL=/bin/bash

local_version := $(shell python3 setup.py --version)
tmpdir := $(shell mktemp -d)

readme: tcdemux/README.rst

tcdemux/README.rst: README.md
	pandoc -f markdown -t rst README.md > tcdemux/README.rst


tests: test_container test_local

test_container:
	$(eval container_version=$(shell wget -qO- https://quay.io/api/v1/repository/biocontainers/tcdemux/tag/ | jq -r '.tags[0].name' ))
	snakemake -s extras/tests.smk \
	--config version=$(container_version) outdir="tests/container_$(container_version)" \
	--use-singularity \
	--singularity-args "-B /datasets,$(PWD),$(tmpdir) -H $(tmpdir) --pwd $(PWD) --containall --cleanenv --writable-tmpfs" \
	--printshellcmds \
	--cores 10 \

test_local:
	snakemake -s extras/tests.smk \
	--config version="" outdir="tests/local_$(local_version)" \
	--cores 10 \

graph: assets/graph.svg assets/external_only_graph.svg

assets/graph.svg:
	snakemake -n \
	-s tcdemux/Snakefile \
	--cores 12 \
	--rulegraph \
	--forceall \
	--config \
	sample_data_file=test-data/samples_graph.csv \
	read_directory=test-data/raw_reads \
	adaptor_files=test-data/adaptors/TruSeq3-PE-2.fa \
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
	sample_data_file=test-data/samples_unpooled_graph.csv \
	read_directory=test-data/raw_reads \
	adaptor_files=test-data/adaptors/TruSeq3-PE-2.fa \
	outdir=tests/graph \
	mem_gb=12 \
	threads=12 \
	keep_intermediate_files=False \
	| grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
	| dot -Tsvg \
	> assets/external_only_graph.svg
