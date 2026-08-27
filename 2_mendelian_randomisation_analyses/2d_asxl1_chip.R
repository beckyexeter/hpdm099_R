library(data.table)
library(ggplot2)
library(ieugwasr)
library(MRPRESSO)
library(showtext)
library(TwoSampleMR)

if (file.exists("results/chip-subtypes/asxl1_chip.txt")) {
    file.remove("results/chip-subtypes/asxl1_chip.txt")
}

font_add_google("Alegreya Sans", family = "alegreya_sans")
showtext_auto()

exp_raw <- fread("data/li_ltl.tsv.gz")
exp_raw <- subset(exp_raw, exp_raw$pval < 5e-8)
exp_dat <- format_data(as.data.frame(exp_raw))
exp_dat$SNP <- sub("^(rs[0-9]+).*", "\\1", exp_dat$SNP)
clumped_exp <- clump_data(exp_dat)
exp_dat <- clumped_exp
exp_dat$exposure <- "telomere length"
rm(exp_raw, clumped_exp)

out_raw <- fread("data/kessler_asxl1_chip.tsv.gz") 
out_raw$SNP <- sub("^(rs[0-9]+).*", "\\1", out_raw$SNP)
out_raw <- subset(out_raw, out_raw$SNP %in% exp_dat$SNP)
out_raw$keep <- FALSE
for (i in 1:nrow(out_raw)) {
    flips <- c(A = "T", T = "A", C = "G", G = "C")
    snp <- out_raw$SNP[i]
    j <- match(snp, exp_dat$SNP)
    a1 <- out_raw$effect_allele[i]
    a2 <- out_raw$other_allele[i]
    b1 <- exp_dat$effect_allele.exposure[j]
    b2 <- exp_dat$other_allele.exposure[j]
    if (a1 == b1 && a2 == b2) {
        out_raw$keep[i] <- TRUE
    } else if (a1 == b2 && a2 == b1) {
        out_raw$keep[i] <- TRUE
    } else if (a1 == flips[b1] && a2 == flips[b2]) {
        out_raw$keep[i] <- TRUE
    } else if (a1 == flips[b2] && a2 == flips[b1]) {
        out_raw$keep[i] <- TRUE
    } else {
        out_raw$keep[i] <- FALSE
    }
}
out_matched <- out_raw[out_raw$keep, ]
out_matched$keep <- NULL 
out_dat <- format_data(as.data.frame(out_matched), type="outcome")
out_dat$outcome <- "ASXL1 CHIP"
rm(out_raw, out_matched, i, flips, snp, j, a1, a2, b1, b2)

harmonised_data <- harmonise_data(exp_dat, out_dat)
harmonised_data <- harmonised_data[harmonised_data$pval.outcome > 5e-8, ]

mr_presso_res <- mr_presso(BetaOutcome = "beta.outcome",
                           BetaExposure = "beta.exposure",
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                           data = harmonised_data, NbDistribution = 1000,
                           SignifThreshold = 0.05)

if (is.na(mr_presso_res$`Main MR results`$`Causal Estimate`[2])) {
    harmonised_data_without_outliers <- harmonised_data
} else {
    harmonised_data_without_outliers <- harmonised_data[-c(mr_presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`), ]
}

instrument_f_statistics <- data.frame(
    snp = harmonised_data_without_outliers$SNP,
    f_statistic = (harmonised_data_without_outliers$beta.exposure / harmonised_data_without_outliers$se.exposure)^2
)

mr_res <- mr(harmonised_data_without_outliers, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")) 

res_pleio <- mr_pleiotropy_test(harmonised_data_without_outliers)

res_het <- mr_heterogeneity(harmonised_data, method_list = c("mr_ivw", "mr_egger_regression"))

method_labels <- c(
    "Inverse variance weighted" = "IVW",
    "MR Egger" = "MR-Egger",
    "Weighted median" = "Weighted median",
    "Weighted mode" = "Weighted mode"
)

pval_labs <- setNames(
    sprintf("%s: P = %.1e", method_labels[mr_res$method], mr_res$pval),
    mr_res$method
)

mr_colours <- c(
    "Inverse variance weighted" = "#a6cee3",
    "MR Egger" = "#1f78b4",
    "Weighted median" = "#b2df8a",
    "Weighted mode" = "#33a02c"
)

p1 <- mr_scatter_plot(mr_res, harmonised_data_without_outliers)[[1]]

p1$layers <- append(
    p1$layers,
    list(
        geom_hline(yintercept = 0, colour = "darkgray"),
        geom_vline(xintercept = 0, colour = "darkgray")
    ),
    after = 0
)

p1 <- p1 + 
    theme_bw() + 
    ylab(expression(paste("SNP effect on ", italic("ASXL1"), "CHIP"))) +
    scale_colour_manual(values = mr_colours, labels = pval_labs) +
    theme(
        text = element_text(family = "alegreya_sans", size = 90),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c("right", "top")
    )

res_single <- mr_singlesnp(harmonised_data_without_outliers, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

p2 <- mr_forest_plot(res_single)[[1]] +
    theme_bw() +
    xlab(expression(paste("MR effect size for telomere length on ", italic("ASXL1"), "CHIP"))) +
    theme(
        text = element_text(family = "alegreya_sans", size = 75),
        legend.position = "none"
    ) 

ggsave("plots/chip-subtypes/asxl1_chip_scatter.png", p1, width = 9, height = 6, dpi = 600)
ggsave("plots/chip-subtypes/asxl1_chip_loo.png", p2, width = 6, height = 6, dpi = 600)

cat("MR Analysis of telomere length on ASXL1 CHIP\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

cat("\nMR-PRESSO Observed Residual Sum of Squares:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(mr_presso_res$`MR-PRESSO results`$`Global Test`$`RSSobs`),
               file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

cat("\nMR-PRESSO Global Test P-value:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(mr_presso_res$`MR-PRESSO results`$`Global Test`$`Pvalue`),
               file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

if (!is.na(mr_presso_res$`Main MR results`$`Causal Estimate`[2])) {
    cat("\nMR-PRESSO Outlier Indices:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
    capture.output(print(mr_presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                   file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
} 

cat("\ntelomere length on CHIP MR Results:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(mr_res), file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

cat("\nMR Egger Pleiotropy Test:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(res_pleio), file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

cat("\nCochran's Q:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(res_het), file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)

cat("\nInstrument F Statistics:\n", file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
capture.output(print(instrument_f_statistics), file = "results/chip-subtypes/asxl1_chip.txt", append = TRUE)
