library(data.table)
library(forestploter)
library(grid)
library(gridtext)
library(ieugwasr)
library(MRPRESSO)
library(showtext)
library(TwoSampleMR)

font_add_google("Alegreya Sans", family = "alegreya_sans")
showtext_auto()

outcomes <- c(
    "CHIP (any subtype)", "DNMT3A CHIP", "TET2 CHIP", "ASXL1 CHIP",
    "PPM1D CHIP", "TP53 CHIP", "JAK2 CHIP", "SF3B1 CHIP", "SRSF2 CHIP"
)

outcome_data_paths <- c("data/kessler_all_chip.tsv.gz", "data/kessler_dnmt3a_chip.tsv.gz",
                        "data/kessler_tet2_chip.tsv.gz", "data/kessler_asxl1_chip.tsv.gz",
                        "data/kessler_ppm1d_chip.tsv.gz", "data/kessler_tp53_chip.tsv.gz",
                        "data/kessler_jak2_chip.tsv.gz", "data/kessler_sf3b1_chip.tsv.gz",
                        "data/kessler_srsf2_chip.tsv.gz")

mr_res <- data.frame(
    outcome = outcomes,
    SNPs = rep(NA,, 9),
    beta = rep(NA,, 9),
    se = rep(NA,, 9),
    p_value = rep(NA,, 9),
    egger_p = rep(NA,, 9),
    mr_presso_rss = rep(NA,, 9),
    mr_presso_global_p = rep(NA,, 9),
    mr_presso_outliers_removed = rep(FALSE,, 9),
    mr_presso_num_outliers_removed = rep(NA,, 9)
)

exp_raw <- fread("data/li_ltl.tsv.gz")
exp_raw <- subset(exp_raw, exp_raw$pval < 5e-8)
exp_dat <- format_data(as.data.frame(exp_raw))
exp_dat$SNP <- sub("^(rs[0-9]+).*", "\\1", exp_dat$SNP)
clumped_exp <- clump_data(exp_dat)
exp_dat <- clumped_exp
exp_dat$exposure <- "telomere length"
rm(exp_raw, clumped_exp)

for (out_num in 1:9) {
    out_raw <- fread(outcome_data_paths[out_num]) 
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
    out_dat$outcome <- outcomes[out_num]
    rm(out_raw, out_matched, i, flips, snp, j, a1, a2, b1, b2)

    harmonised_data <- harmonise_data(exp_dat, out_dat)
    harmonised_data <- harmonised_data[harmonised_data$pval.outcome > 5e-8, ]

    mr_presso_res <- mr_presso(BetaOutcome = "beta.outcome",
                               BetaExposure = "beta.exposure",
                               SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                               data = harmonised_data, NbDistribution = 1000,
                               SignifThreshold = 0.05)

    mr_res$mr_presso_rss[out_num] <- mr_presso_res$`MR-PRESSO results`$`Global Test`$`RSSobs`

    mr_res$mr_presso_global_p[out_num] <- mr_presso_res$`MR-PRESSO results`$`Global Test`$`Pvalue`

    if (is.na(mr_presso_res$`Main MR results`$`Causal Estimate`[2])) {
        harmonised_data_without_outliers <- harmonised_data
    } else {
        harmonised_data_without_outliers <- harmonised_data[-c(mr_presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`), ]
        mr_res$mr_presso_outliers_removed[out_num] <- TRUE
        mr_res$mr_presso_num_outliers_removed[out_num] <- length(mr_presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) 
    }

    mr_out <- mr(harmonised_data_without_outliers, method_list = c("mr_ivw", "mr_egger_regression"))
    res_pleio <- mr_pleiotropy_test(harmonised_data_without_outliers)

    mr_res$SNPs[out_num] <- mr_out$nsnp[1]
    mr_res$beta[out_num] <- mr_out$b[1]
    mr_res$se[out_num] <- mr_out$se[1]
    mr_res$p_value[out_num] <- mr_out$pval[1]
    mr_res$egger_p[out_num] <- res_pleio$pval[1]
}

odds_ratios <- exp(mr_res$beta)
lowers <- exp(mr_res$beta - qnorm(0.975) * mr_res$se)
uppers <- exp(mr_res$beta + qnorm(0.975) * mr_res$se)

plot_tbl <- data.frame(
    Outcome = paste(rep(" ", 50), collapse = ""),
    SNPs = mr_res$SNPs,
    `OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", odds_ratios, lowers, uppers),
    `    ` = paste(rep(" ", 100), collapse = " "),
    `P-value` = ifelse(mr_res$p_value < 0.001, "<0.001", sprintf("%.3f", mr_res$p_value)),
    `Egger intercept\nP-value` = sprintf("%.3f", mr_res$egger_p),
    check.names = FALSE
)

p <- forest(
    plot_tbl,
    est = odds_ratios,
    lower = lowers,
    upper = uppers,
    ci_column = 4,
    ref_line = 1,
    xlim = c(0.1, 10),
    ticks_at = c(0.1, 1, 10),
    x_trans = "log",
    theme = forest_theme(
        base_family = "alegreya_sans",
        base_size = 9,
        ci_pch = 15,
        ci_col = "black",
        refline_lty = "dashed",
        refline_col = "grey40"
    )
)

p <- edit_plot(p, col = c(5,6), 
               which = "text",
               hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

gene_names <- c(NA, "DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "JAK2", "SF3B1", "SRSF2")

for (i in 1:9) {
    colour <- if (i <= 3) "blue" else "black"
    
    label <- if (i == 1) {
        "CHIP (any subtype)"
    } else {
        bquote(italic(.(gene_names[i])) ~ CHIP)
    }
    
    p <- add_grob(
        p,
        row = i, col = 1,
        order = "top",
        gb_fn = textGrob,
        label = label,
        x = unit(4, "pt"), y = unit(0.5, "npc"),
        hjust = 0, vjust = 0.5,
        gp = gpar(fontsize = 9, fontfamily = "alegreya_sans", col = colour)
    )
}

dims <- get_wh(plot = p, unit = "in")

png(
    filename = "plots/ltl_chip_mr_forest_plot.png",
    width = dims[1],
    height = dims[2],
    units = "in",
    res = 600
)

par(mar = c(0, 0, 0, 0))
grid.newpage()
grid.draw(p)
dev.off()

write.table(mr_res, "results/ltl_chip_mr_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")
