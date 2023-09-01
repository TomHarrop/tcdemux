readme: README.rst

README.rst: README.md
	pandoc -f markdown -t rst README.md > README.rst

test: pool_test/134567_LibID134667_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S7_L002_r1.fastq.gz nopool_test/134567_LibID134667_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S7_L002_r1.fastq.gz

pool_test/134567_LibID134667_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S7_L002_r1.fastq.gz:
	tcdemux \
	--sample_data data/samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir pool_test \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \

nopool_test/134567_LibID134667_GAP_BRF_H5TT7DRX3_CCGCGGTT-CTAGCGCT_S7_L002_r1.fastq.gz:
	tcdemux \
	--sample_data data/unpooled_samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir nopool_test \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16

insane:
	tcdemux \
	--sample_data data/insane_samples.csv \
	--adaptors data/adaptors/alicia_adapters.fa data/adaptors/TruSeq3-PE-2.fa data/adaptors/bbmap_39.01_adaptors.fa \
	--outdir insane_test \
	--read_directory data/test_data \
	--threads 10 \
	--mem_gb 16 \

