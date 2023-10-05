#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd

stats_file = Path("tests/pooled/stats/cutadapt/AG_P8P5C.json")

with stats_file.open("r") as f:
    stats_json = json.loads(f.read())

# create dataframe
cutadapt_stats = pd.json_normalize(stats_json)

# parse the adaptor info
library_count = {}
adaptor_stats = cutadapt_stats.iloc[0]["adapters_read1"]

for sample in adaptor_stats:
    library_count[sample["name"]] = sample["total_matches"]

s = pd.Series(library_count, name="reads_assigned")
s.index.name = "library"
processed_stats = pd.DataFrame(s)
processed_stats["total_input"] = cutadapt_stats.iloc[0]["read_counts.input"]
processed_stats["total_matched"] = cutadapt_stats.iloc[0][
    "read_counts.read1_with_adapter"
]

processed_stats["fraction_of_total"] = processed_stats["reads_assigned"].div(
    processed_stats["total_input"]
)

processed_stats["fraction_of_matched"] = processed_stats["reads_assigned"].div(
    processed_stats["total_matched"]
)
