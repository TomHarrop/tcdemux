#!/usr/bin/env Rscript

library(readxl)
library(data.table)

metadata_file <- "data/config_file_prep/2022-03-03_GAP_TEST_METADATA.xlsx"
raw_data_dir <- "data/raw_data"

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

all_pools <- all_metadata[, unique(pool_name)]
all_pools <- as.character(na.omit(all_pools))
names(all_pools) <- all_pools

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

pool_files <- rbindlist(
  lapply(all_pools, GetReadFiles),
  idcol = "pool_name")
