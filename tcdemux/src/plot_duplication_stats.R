#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)
library(ggplot2)

before_file <- snakemake@input[["before_file"]]
after_file <- snakemake@input[["after_file"]]
plot_file <- snakemake@output[["plot"]]

dup_data_long <- rbindlist(
  list(
    Raw = fread(before_file),
    Processed = fread(after_file)
  ),
  idcol = "step")

dup_data <- dcast(
  dup_data_long,
  sample + step ~ type,
  value.var = "reads")

dup_data[, "Duplication rate" := Duplicates / Input]
dup_data[, step := factor(step, levels = c("Raw", "Processed"))]

gp <- ggplot(dup_data, aes(x = `Duplication rate`)) +
  theme_grey(base_size = 8) +
  facet_grid(. ~ step, scales = "fixed") +
  ylab("Number of samples") +
  geom_histogram(binwidth = 0.01)

ggsave(plot_file,
       gp,
       width = 13.33,
       height = 7.5,
       units = "in",
       device = cairo_pdf,
       limitsize = FALSE)

sessionInfo()

