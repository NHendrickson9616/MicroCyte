library(ggplot2)
library(viridis)
source("bin/scripts/analysis/icellate.R")
source("bin/scripts/analysis/SirMixaPlot.R")

ergodic_slicing <- function(df = cells,
                            subset_df = cells[cells$edu == "Positive", ],
                            variable = "dna_norm",
                            bins = 10,
                            bin_add_name = "",
                            rate = 24,
                            rate_description = "doubling time (hours)",
                            icell = FALSE,
                            icell_number = 3,
                            icell_fusion = FALSE,
                            pca = FALSE,
                            pca_groups = FALSE,
                            hist_save = "figures/ergodicHist.pdf",
                            ergo_save = "figures/ergodicRates.pdf",
                            ergo_file_save = "data/ergo.csv") {
  list_df <- split(subset_df, cut(unlist(subset_df[variable]), breaks = bins))
  if (icell == TRUE) {
    tick <- 1
    for (bin_set in list_df) {
      if (nrow(bin_set) > 0) {
        if (nrow(bin_set) < icell_number) {
          johnny <- nrow(bin_set)
        } else {
          johnny <- icell_number
        }
        icellate(
          target_cells = bin_set,
          folder_name = paste0(bin_add_name, "bin ", tick),
          verify_size = FALSE,
          fuse = icell_fusion,
          randomize = TRUE,
          sampling_number = johnny,
          line_analyses = FALSE
        )
        tick <- tick + 1
      }
    }
  }

  histo_draft <- ggplot(data = subset_df, aes(x = unlist(subset_df[variable]))) +
    geom_density() +
    xlab(variable) +
    theme_classic()
  print(histo_draft)
  ggsave(hist_save)

  ergo_df <- data.frame(
    "bin_id" = c(1:bins),
    "bin_rate" = 0,
    "volume" = 0,
    "mean_stat" = 0
  )
  tick <- 1
  for (bin_set in list_df) {
    if (pca == TRUE) {
      if (!"pcas" %in% list.files(path = "figures/")) {
        dir.create("figures/pcas")
      }
      runPCA(
        df = bin_set,
        groupz = pca_groups,
        save_file = paste0("figures/pcas/", bin_add_name, "bin_", tick, ".png")
      )
    }
    ergo_df[tick, ]$volume <- nrow(bin_set)
    ergo_df[tick, ]$bin_rate <- rate * (
      (2 - (nrow(subset_df) / nrow(df))) / (nrow(bin_set) / nrow(subset_df))
    )
    ergo_df[tick, ]$mean_stat <- median(unlist(bin_set[variable]))
    tick <- tick + 1
  }
  ergo_df$bin_rate <- ergo_df$bin_rate / mean(ergo_df[!is.infinite(ergo_df$bin_rate), ]$bin_rate)
  cbins <- c(1:bins)
  if (TRUE %in% is.na(ergo_df$mean_stat)) {
    ergo_df[is.na(ergo_df$mean_stat), ]$mean_stat <- ""
  }
  ergo_draft <- ggplot(
    data = ergo_df[!is.infinite(ergo_df$bin_rate), ],
    aes(x = bin_id, y = bin_rate, fill = volume)
  ) +
    geom_bar(stat = "identity", position = position_dodge(0.7)) +
    ylab("Mean bin progression rate (a.u. per hour)") +
    #  scale_x_continuous(breaks = cbins, labels = ergo_df$mean_stat)+
    xlab(paste0("Median size of ", variable)) +
    theme_classic() +
    scale_fill_viridis() +
    theme(legend.position = "top")
  print(ergo_draft)
  ggsave(ergo_save)
  write.csv(ergo_df, file = ergo_file_save, row.names = FALSE)
}
