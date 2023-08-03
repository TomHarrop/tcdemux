#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

sample_data_file = snakemake.input['sample_data']
barcode_file = snakemake.output['barcode_file']
mypool = snakemake.params['pool']

sample_data = pd.read_csv(
    sample_data_file, index_col='Name')

# subset the data
pool_sd = sample_data[sample_data['pool_name'] == mypool]
pool_samples = sorted(set(pool_sd.index))
# prepare the adaptor sequences
seq_records = []
for idx, row in pool_sd.iterrows():
    seq = Seq(row['internal_index_sequence'])
    record = SeqRecord(
        seq,
        id=str(idx),
        description=''
        )
    seq_records.append(record)
# write the records to the barcode file
with open(barcode_file, 'w') as f:
    SeqIO.write(seq_records, f, 'fasta')
