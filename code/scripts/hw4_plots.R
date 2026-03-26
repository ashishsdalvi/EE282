library(ggplot2)

in_dir <- "/pub/asdalvi/informatics_class/ee282_hw4/data/processed/"
out_dir <- "/pub/asdalvi/informatics_class/ee282_hw4/output/figures/"

d_large <- read.table(paste0(in_dir, "data_large.txt"), col.names=c("len", "gc"))
d_small <- read.table(paste0(in_dir, "data_smalll.txt"), col.names=c("len", "gc"))

make_plots <- function(df, suffix) {
  p1 <- ggplot(df, aes(x=len)) + geom_histogram(bins=50) + scale_x_log10() + theme_bw()
  ggsave(paste0(out_dir, "hist_len_", suffix, ".png"), p1)
  
  p2 <- ggplot(df, aes(x=gc)) + geom_histogram(bins=50) + theme_bw()
  ggsave(paste0(out_dir, "hist_gc_", suffix, ".png"), p2)
  
  df_sorted <- df[order(-df$len), ]
  df_sorted$cum_size <- cumsum(as.numeric(df_sorted$len))
  df_sorted$rank <- 1:nrow(df_sorted)
  p3 <- ggplot(df_sorted, aes(x=rank, y=cum_size)) + geom_line() + theme_bw()
  ggsave(paste0(out_dir, "cdf_", suffix, ".png"), p3)
}

make_plots(d_large, "large")
make_plots(d_small, "smalll")
