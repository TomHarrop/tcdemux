#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)
library(ggplot2)

trim_files <- snakemake@input[["trim_files"]]
summary_file <- snakemake@input[["summary_file"]]
plot_file <- snakemake@output[["plot"]]


summary_data <- fread(summary_file)

names(trim_files) <- sapply(trim_files, function(x) 
  unlist(strsplit(basename(x), ".", fixed = TRUE))[[1]])

trim_list <- lapply(trim_files, function(x)
  fread(x, fill = TRUE, skip = 3)[, .(adaptor = `#Name`, reads = Reads)])
trim_data <- rbindlist(trim_list, idcol = "sample")

# get the total reads per indiv
summary_wide <- dcast(summary_data,
                      sample ~ type,
                      value.var = c("reads", "bases"))

trim_with_summary <- merge(trim_data,
                           summary_wide,
                           by = "sample")

trim_with_summary[, adaptor_frac := reads / reads_Input ]

# get the top5 adaptors for each indiv
setorder(trim_with_summary, -adaptor_frac)
top5i <- trim_with_summary[, .I[1:5], by = sample][, V1]
adaptors_to_plot <- trim_with_summary[top5i, unique(adaptor)]
adaptor_order <- trim_with_summary[, sum(reads), by = adaptor][order(-V1), unique(adaptor)]

# add the fraction of bases trimmed
trim_with_summary[, bases_result_frac := bases_Result / bases_Input]

# subset the plot data
adaptor_pd <- trim_with_summary[adaptor %in% adaptors_to_plot]
adaptor_pd[, adaptor := factor(adaptor, levels = rev(adaptor_order))]

# go long for facet_plot
adaptor_pd_long <- melt(adaptor_pd,
                        id.vars = c("sample", "adaptor"),
                        measure.vars = c("adaptor_frac", "bases_result_frac", "reads_Result"))
adaptor_pd_long[variable == "adaptor_frac",
                facet_lab := "Fraction of reads with adaptor"]
adaptor_pd_long[variable == "bases_result_frac",
                facet_lab := "Fraction of bases after trimming"]
adaptor_pd_long[variable == "reads_Result",
                facet_lab := "Reads after trimming (M)"]

adaptor_pd_long[, facet_lab := relevel(factor(facet_lab),
                                       "Fraction of reads with adaptor")]

# get order for indivs
sample_order <- trim_with_summary[order(-reads_Result), unique(sample)]
adaptor_pd_long[, sample := factor(sample, levels = sample_order)]

# plot with both facets
extracols <- viridis::viridis(3, option = "A")
gp <- ggplot(mapping = aes(x = sample, y = value, fill = adaptor)) +
  theme_grey(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = "top",
        legend.key.size = unit(8, "pt")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(guide = guide_legend(title = NULL,
                                            reverse = TRUE)) +
  xlab(NULL) + ylab(NULL) +
  facet_grid(. ~ facet_lab,
             scales = "free_x",
             switch = "x") +
  coord_flip() +
  geom_col(data = adaptor_pd_long[variable == "adaptor_frac"]) + 
  geom_col(data = unique(adaptor_pd_long[variable == "bases_result_frac"],
                         by = "sample"),
           fill = extracols[[1]]) +
  geom_col(data = unique(adaptor_pd_long[variable == "reads_Result"],
                         by = "sample"),
           aes(y = value / 1e6),
           fill = extracols[[2]])

# 3mm per line plus 30 mm
ho_mm <- (length(sample_order) * 3) + 30

ggsave(plot_file,
       gp,
       width = 210,
       height = ho_mm,
       units = "mm",
       device = cairo_pdf,
       limitsize = FALSE)

sessionInfo()
