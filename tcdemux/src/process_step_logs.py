#!/usr/bin/env python3

import sys
import pandas as pd

# Check if the sample name argument is provided
if len(sys.argv) != 2:
    print("Usage: script.py sample_name", file=sys.stderr)
    sys.exit(1)

sample_name = sys.argv[1]

# read from stdin
df = pd.read_csv(sys.stdin, delimiter="\t", header=None)

# mung the columns
df.columns = ["V1", "V2", "V3"]
df["type"] = df["V1"].str.replace(r"[^a-zA-Z]", "", regex=True)
df["reads"] = pd.to_numeric(df["V2"].str.split().str[0])
df["bases"] = pd.to_numeric(df["V3"].str.split().str[0])

# Add the sample name as an extra column
df["sample"] = sample_name

# write to stdout
df[["sample", "type", "reads", "bases"]].to_csv(sys.stdout, index=False, header=False)
