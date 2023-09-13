readme: README.rst

README.rst: README.md
	pandoc -f markdown -t rst README.md > README.rst

test: pooled unpooled dirty missing

pooled: tests/pooled/134501_LibID134601_GAP_BRF_H5TT7DRX3_GCTTGTCA-GTATGTTC_S82_L002_r1.fastq.gz
tests/pooled/134501_LibID134601_GAP_BRF_H5TT7DRX3_GCTTGTCA-GTATGTTC_S82_L002_r1.fastq.gz:
	tcdemux \
	--sample_data data/samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/pooled \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \

unpooled: tests/unpooled/134486_LibID134586_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S10_L002_r1.fastq.gz
tests/unpooled/134486_LibID134586_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S10_L002_r1.fastq.gz:
	tcdemux \
	--sample_data data/unpooled_samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/unpooled \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16

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
	tcdemux \
	--sample_data data/samples_missing.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir tests/missing \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \

