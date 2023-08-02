#!/usr/bin/env Rscript

library(readxl)
library(data.table)

metadata_file <- "data/config_file_prep/2022-03-03_GAP_TEST_METADATA.xlsx"
raw_data_dir <- "data/raw_data"
processed_config_file <- "data/samples.csv"

all_metadata <- as.data.table(
  readxl::read_xlsx(
    metadata_file, sheet = "CSIRO_method"
  )
)

# find the pool name fields
dput(names(all_metadata)[startsWith(names(all_metadata), "Pool")])
pool_name_fields <- c("Pool name...20", "Pool name...64")

# checks
all_metadata[, ..pool_name_fields]
all_metadata[`Pool name...20` == `Pool name...64`]

# just use the first field that matches
pool_name_field <- names(
  all_metadata
)[startsWith(
  names(all_metadata),
  "Pool name")][[1]]

all_metadata[, "pool_name" := get(pool_name_field)]

# find out which barcode combo should be in each pool
i5_colname <- names(all_metadata)[
  startsWith(
    names(all_metadata),
    "i5 external index sequence"
  )][[1]]
i7_colname <- names(all_metadata)[
  startsWith(
    names(all_metadata),
    "i7 external index sequence"
  )][[1]]

all_metadata[, "i5_index" := get(i5_colname)]
all_metadata[, "i7_index" := get(i7_colname)]

# subset the data for external barcode demultiplexing
index_data <- unique(
  all_metadata[, .(pool_name, i5_index, i7_index)]
)
# hopefully, NAs are just empty rows from excel
index_data <- na.omit(index_data)

# try to find the pools in the read directory
GetReadFiles <- function(pool_name){
  all_read_files <- list.files(
    raw_data_dir,
    pattern = pool_name
  )
  r1_file <- grep("R1.*f(ast)?q",
                  all_read_files,
                  value = TRUE,
                  ignore.case = TRUE)
  r2_file <- grep("R2.*f(ast)?q",
                  all_read_files,
                  value = TRUE,
                  ignore.case = TRUE)
  return(
    data.table(
      r1_file = r1_file,
      r2_file = r2_file
    )
  )
}

index_data[
  ,
  c("r1_file", "r2_file") := GetReadFiles(pool_name),
  by = pool_name]

# clean up some colnames
all_metadata[, internal_index_sequence := `Internal index sequence`]

fn_field <- names(all_metadata)[
  startsWith(
    names(all_metadata),
    "Filename"
  )][[1]]
all_metadata[, "filename" := get(fn_field)]

# now get the sample data, including the pool name
cols <- c("Name",
          "filename",
          "internal_index_sequence",
          "pool_name",
          "i5_index",
          "i7_index"
          )

tidy_metadata <- na.omit(all_metadata[, ..cols])

# clean up the sample names (TODO: get human readable sample names)
tidy_metadata[, Name := sub("_bbduk", "", Name)]

config_file <- merge(
  tidy_metadata,
  index_data,
  by = c("pool_name", "i5_index", "i7_index"),
  all.x = TRUE,
  all.y = FALSE)

fwrite(config_file, processed_config_file)

sessionInfo()
