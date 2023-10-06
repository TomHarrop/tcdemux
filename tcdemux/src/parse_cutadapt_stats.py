#!/usr/bin/env python3

import sys
import json
import pandas as pd

# Read JSON input from stdin
input_json = sys.stdin.read()
stats_json = json.loads(input_json)

# Create a DataFrame from the JSON data
cutadapt_stats = pd.json_normalize(stats_json)

# Extract adapter statistics for each library
library_count = {}
adaptor_stats = cutadapt_stats.iloc[0]["adapters_read1"]

for sample in adaptor_stats:
    library_count[sample["name"]] = sample["total_matches"]

# Create a Series with library counts and assign it to a DataFrame
processed_stats = pd.Series(library_count, name="reads_assigned").reset_index()
processed_stats.rename(columns={"index": "library"}, inplace=True)

# Add additional columns to the DataFrame
processed_stats["pool_total_input"] = cutadapt_stats.iloc[0]["read_counts.input"]
processed_stats["pool_total_matched"] = cutadapt_stats.iloc[0][
    "read_counts.read1_with_adapter"
]
processed_stats["fraction_of_total"] = (
    processed_stats["reads_assigned"] / processed_stats["pool_total_input"]
)
processed_stats["fraction_of_matched"] = (
    processed_stats["reads_assigned"] / processed_stats["pool_total_matched"]
)

# Output the processed CSV data to stdout
processed_stats.to_csv(sys.stdout, index=False)
