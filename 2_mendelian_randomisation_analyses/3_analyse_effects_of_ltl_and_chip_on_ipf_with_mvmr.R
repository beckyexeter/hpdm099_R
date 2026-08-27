library(data.table)
library(ggplot2)
library(glmnet)
library(ieugwasr)
library(MASS)
library(MendelianRandomization)
library(MVMR)
library(MVMRmode)
library(MRPRESSO)
library(patchwork)
library(quantreg)
library(robustbase)
library(showtext)
library(TwoSampleMR)

update_geom_defaults("point", list(shape = 4, colour = "blue", size = 2))
font_add_google("Alegreya Sans", family = "alegreya_sans")
showtext_auto()

results_estimates <- data.frame(
    outcome_data = c(rep("gbmi", 4), rep("FinnGen", 4)),
    method = rep(c("MRMR IVW", "MVMR Egger", "MVMR WMe", "MVMR WMo"), 2),
    ltl_beta = rep(NA, 8),
    ltl_se = rep(NA, 8),
    ltl_p_value = rep(NA, 8),
    ltl_odds_ratio = rep(NA, 8),
    ltl_or_lower = rep(NA, 8),
    ltl_or_upper = rep(NA, 8),
    chip_beta = rep(NA, 8),
    chip_se = rep(NA, 8),
    chip_p_value = rep(NA, 8),
    chip_odds_ratio = rep(NA, 8),
    chip_or_lower = rep(NA, 8),
    chip_or_upper = rep(NA, 8)
)

results_pleiotropy <- data.frame(
    outcome_data = c("GBMI", "FinnGen"),
    heter_stat = rep(NA, 2),
    heter_p = rep(NA, 2),
    mr_presso_rss = rep(NA, 2),
    mr_presso_global_test_p = rep(NA, 2),
    mr_presso_num_outliers = rep(NA, 2),
    egger_intercept = rep(NA, 2),
    egger_intercept_p = rep(NA, 2)
)

results_instruments <- data.frame(
    outcome_data = c("GBMI", "FinnGen"),
    total_snps = rep(NA, 2),
    ltl_snps = rep(NA, 2),
    chip_snps = rep(NA, 2),
    overlap = rep(NA, 2)
)

results_instrument_strength <- data.frame(
    outcome_data = c("GBMI", "FinnGen"),
    ltl_cond_f = rep(NA, 2),
    chip_cond_f = rep(NA, 2)
)

chip_exp <- fread("data/kessler_all_chip.tsv.gz")
chip_exp_filtered <- subset(chip_exp, chip_exp$pval < 5e-8)
chip_exp_filtered <- format_data(as.data.frame(chip_exp_filtered))
chip_exp_filtered$SNP <- sub("^(rs[0-9]+).*", "\\1", chip_exp_filtered$SNP)
chip_instruments <- clump_data(chip_exp_filtered)
rm(chip_exp_filtered)

ltl_exp <- fread("data/li_ltl.tsv.gz")
ltl_exp_filtered <- subset(ltl_exp, ltl_exp$pval < 5e-8)
ltl_exp_filtered <- format_data(as.data.frame(ltl_exp_filtered))
ltl_exp_filtered$SNP <- sub("^(rs[0-9]+).*", "\\1", ltl_exp_filtered$SNP)
ltl_instruments <- clump_data(ltl_exp_filtered)
rm(ltl_exp_filtered)

chip_exp$SNP <- sub("^(rs[0-9]+).*", "\\1", chip_exp$SNP)
chip_data_for_ltl_instruments <- subset(chip_exp, chip_exp$SNP %in% ltl_instruments$SNP)
chip_data_for_ltl_instruments <- format_data(as.data.frame(chip_data_for_ltl_instruments))

ltl_exp$SNP <- sub("^(rs[0-9]+).*", "\\1", ltl_exp$SNP)
ltl_data_for_chip_instruments <- subset(ltl_exp, ltl_exp$SNP %in% chip_instruments$SNP)
ltl_data_for_chip_instruments <- format_data(as.data.frame(ltl_data_for_chip_instruments))

chip_exp <- rbind(chip_instruments, chip_data_for_ltl_instruments)
chip_exp <- format_data(
    chip_exp,
    type = "exposure",
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    chr_col = "chr.exposure"
)
chip_exp$exposure <- "CHIP"

ltl_exp <- rbind(ltl_instruments, ltl_data_for_chip_instruments)
ltl_exp <- format_data(
    ltl_exp,
    type = "exposure",
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    chr_col = "chr.exposure"
)
ltl_exp$exposure <- "telomere length"

comb_exp <- rbind(chip_exp, ltl_exp)
comb_exp$exposure <- "EXPOSURE"
comb_exp <- format_data(
    comb_exp,
    type = "exposure",
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    chr_col = "chr.exposure"
)
comb_exp$exposure <- "EXPOSURE"
clump_comb_exp <- clump_data(comb_exp)

comb_exp <- rbind(chip_exp[chip_exp$SNP %in% clump_comb_exp$SNP, ], ltl_exp[ltl_exp$SNP %in% clump_comb_exp$SNP, ])

gbmi_out_raw <- fread("data/gbmi_eur_ipf.tsv.gz") 
gbmi_out_raw <- subset(gbmi_out_raw, gbmi_out_raw$SNP %in% comb_exp$SNP)
gbmi_out_dat <- format_data(as.data.frame(gbmi_out_raw), type="outcome")
gbmi_out_dat$outcome <- "IPF (GBMI meta-analysis)"
gbmi_out_dat <- gbmi_out_dat[gbmi_out_dat$pval.outcome > 5e-8, ]

fg_out_raw <- fread("data/finngen_ipf.tsv.gz") 
fg_out_raw <- subset(fg_out_raw, fg_out_raw$SNP %in% comb_exp$SNP)
fg_out_dat <- format_data(as.data.frame(fg_out_raw), type="outcome")
fg_out_dat$outcome <- "IPF (FinnGen GWAS)"
fg_out_dat <- fg_out_dat[fg_out_dat$pval.outcome > 5e-8, ]

mvdat_gbmi <- mv_harmonise_data(comb_exp, gbmi_out_dat)

mvmr_dat_gbmi <- cbind(mvdat_gbmi$exposure_beta, 
                       mvdat_gbmi$exposure_se, 
                       mvdat_gbmi$outcome_beta, 
                       mvdat_gbmi$outcome_se
)
names_gbmi <- rownames(mvmr_dat_gbmi)
rownames(mvmr_dat_gbmi) <- NULL
mvmr_dat_gbmi <- cbind(names_gbmi, mvmr_dat_gbmi)
mvmr_dat_gbmi <- as.matrix(mvmr_dat_gbmi)

F_data_gbmi <- format_mvmr(BXGs = mvmr_dat_gbmi[, c(2,3)],
                           BYG = mvmr_dat_gbmi[,6],
                           seBXGs = mvmr_dat_gbmi[,c(4,5)],
                           seBYG = mvmr_dat_gbmi[,7],
                           RSID = mvmr_dat_gbmi[,1])

bx_matrix_gbmi <- as.matrix(cbind(F_data_gbmi$betaX1, F_data_gbmi$betaX2))

sebx_matrix_gbmi <- as.matrix(cbind(bxse = F_data_gbmi$sebetaX1, bxse = F_data_gbmi$sebetaX2))

mvinput_gbmi <- mr_mvinput(bx = bx_matrix_gbmi, bxse = sebx_matrix_gbmi,
                           by = F_data_gbmi$betaYG, byse = F_data_gbmi$sebetaYG)

mvdat_fg <- mv_harmonise_data(comb_exp, fg_out_dat)

mvmr_dat_fg <- cbind(mvdat_fg$exposure_beta, 
                     mvdat_fg$exposure_se, 
                     mvdat_fg$outcome_beta, 
                     mvdat_fg$outcome_se
)
names_fg <- rownames(mvmr_dat_fg)
rownames(mvmr_dat_fg) <- NULL
mvmr_dat_fg <- cbind(names_fg, mvmr_dat_fg)
mvmr_dat_fg <- as.matrix(mvmr_dat_fg)

F_data_fg <- format_mvmr(BXGs = mvmr_dat_fg[,c(2,3)],
                                            BYG = mvmr_dat_fg[,6],
                                            seBXGs = mvmr_dat_fg[,c(4,5)],
                                            seBYG = mvmr_dat_fg[,7],
                                            RSID = mvmr_dat_fg[,1])

bx_matrix_fg <- as.matrix(cbind(F_data_fg$betaX1, F_data_fg$betaX2))

sebx_matrix_fg <- as.matrix(cbind(bxse = F_data_fg$sebetaX1, bxse = F_data_fg$sebetaX2))

mvinput_fg <- mr_mvinput(bx = bx_matrix_fg, bxse = sebx_matrix_fg,
                         by = F_data_fg$betaYG, byse = F_data_fg$sebetaYG)

mr_presso_res_gbmi <- mr_presso(BetaOutcome = "betaYG", BetaExposure = c("betaX1", "betaX2"),
                                SdOutcome = "sebetaYG", SdExposure = c("sebetaX1", "sebetaX2"),
                                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = F_data_gbmi,
                                NbDistribution = 1000, SignifThreshold = 0.05)

mvmr_pleio_res_gbmi <- pleiotropy_mvmr(r_input = F_data_gbmi, gencov = 0)
sre_gbmi <- strength_mvmr(r_input = F_data_gbmi, gencov = 0)

mvmr_ivw_gbmi <- mr_mvivw(mvinput_gbmi)
mvmr_egger_gbmi <- mr_mvegger(mvinput_gbmi)
mvmr_wme_gbmi <- mr_mvmedian(mvinput_gbmi)
mvmr_wmo_gbmi <- mv_mrmode(Bout = mvdat_gbmi$outcome_beta, 
                           Bexp = mvdat_gbmi$exposure_beta,
                           SEout = mvdat_gbmi$outcome_se, 
                           SEexp = mvdat_gbmi$exposure_se, 
                           Psi=0, CIMin = NA, CIMax = NA, CIStep = 0.001, alpha = 0.05,
                           residual="IVW", Mode="CM", weighting = "weighted", stderror = "simple", 
                           phi = 1, distribution = "normal", iterations = 10000)

ltl_instruments_filtered_gbmi <- ltl_instruments[ltl_instruments$SNP %in% gbmi_out_dat$SNP, ]
chip_instruments_filtered_gbmi <- chip_instruments[chip_instruments$SNP %in% gbmi_out_dat$SNP, ]
instruments_intersection_gbmi <- intersect(ltl_instruments_filtered_gbmi$SNP, chip_instruments_filtered_gbmi$SNP)

results_estimates$ltl_beta[1] <- mvmr_ivw_gbmi$Estimate[2]
results_estimates$ltl_se[1] <- mvmr_ivw_gbmi$StdError[2]
results_estimates$ltl_p_value[1] <- mvmr_ivw_gbmi$Pvalue[2]
results_estimates$ltl_odds_ratio[1] <- exp(mvmr_ivw_gbmi$Estimate[2])
results_estimates$ltl_or_lower[1] <- exp(mvmr_ivw_gbmi$Estimate[2] - qnorm(0.975) * mvmr_ivw_gbmi$StdError[2])
results_estimates$ltl_or_upper[1] <- exp(mvmr_ivw_gbmi$Estimate[2] + qnorm(0.975) * mvmr_ivw_gbmi$StdError[2])

results_estimates$chip_beta[1] <- mvmr_ivw_gbmi$Estimate[1]
results_estimates$chip_se[1] <- mvmr_ivw_gbmi$StdError[1]
results_estimates$chip_p_value[1] <- mvmr_ivw_gbmi$Pvalue[1]
results_estimates$chip_odds_ratio[1] <- exp(mvmr_ivw_gbmi$Estimate[1] * log(2))
results_estimates$chip_or_lower[1] <- exp(mvmr_ivw_gbmi$Estimate[1] * log(2) - qnorm(0.975) * mvmr_ivw_gbmi$StdError[1] * log(2))
results_estimates$chip_or_upper[1] <- exp(mvmr_ivw_gbmi$Estimate[1] * log(2) + qnorm(0.975) * mvmr_ivw_gbmi$StdError[1] * log(2))

results_estimates$ltl_beta[2] <- mvmr_egger_gbmi$Estimate[2]
results_estimates$ltl_se[2] <- mvmr_egger_gbmi$StdError.Est[2]
results_estimates$ltl_p_value[2] <- mvmr_egger_gbmi$Pvalue.Est[2]
results_estimates$ltl_odds_ratio[2] <- exp(mvmr_egger_gbmi$Estimate[2])
results_estimates$ltl_or_lower[2] <- exp(mvmr_egger_gbmi$Estimate[2] - qnorm(0.975) * mvmr_egger_gbmi$StdError.Est[2])
results_estimates$ltl_or_upper[2] <- exp(mvmr_egger_gbmi$Estimate[2] + qnorm(0.975) * mvmr_egger_gbmi$StdError.Est[2])

results_estimates$chip_beta[2] <- mvmr_egger_gbmi$Estimate[1]
results_estimates$chip_se[2] <- mvmr_egger_gbmi$StdError.Est[1]
results_estimates$chip_p_value[2] <- mvmr_egger_gbmi$Pvalue.Est[1]
results_estimates$chip_odds_ratio[2] <- exp(mvmr_egger_gbmi$Estimate[1] * log(2))
results_estimates$chip_or_lower[2] <- exp(mvmr_egger_gbmi$Estimate[1] * log(2) - qnorm(0.975) * mvmr_egger_gbmi$StdError.Est[1] * log(2))
results_estimates$chip_or_upper[2] <- exp(mvmr_egger_gbmi$Estimate[1] * log(2) + qnorm(0.975) * mvmr_egger_gbmi$StdError.Est[1] * log(2))

results_estimates$ltl_beta[3] <- mvmr_wme_gbmi$Estimate[2]
results_estimates$ltl_se[3] <- mvmr_wme_gbmi$StdError[2]
results_estimates$ltl_p_value[3] <- mvmr_wme_gbmi$Pvalue[2]
results_estimates$ltl_odds_ratio[3] <- exp(mvmr_wme_gbmi$Estimate[2])
results_estimates$ltl_or_lower[3] <- exp(mvmr_wme_gbmi$Estimate[2] - qnorm(0.975) * mvmr_ivw_gbmi$StdError[2])
results_estimates$ltl_or_upper[3] <- exp(mvmr_wme_gbmi$Estimate[2] + qnorm(0.975) * mvmr_ivw_gbmi$StdError[2])

results_estimates$chip_beta[3] <- mvmr_wme_gbmi$Estimate[1]
results_estimates$chip_se[3] <- mvmr_wme_gbmi$StdError[1]
results_estimates$chip_p_value[3] <- mvmr_wme_gbmi$Pvalue[1]
results_estimates$chip_odds_ratio[3] <- exp(mvmr_wme_gbmi$Estimate[1] * log(2))
results_estimates$chip_or_lower[3] <- exp(mvmr_wme_gbmi$Estimate[1] * log(2) - qnorm(0.975) * mvmr_wme_gbmi$StdError[1] * log(2))
results_estimates$chip_or_upper[3] <- exp(mvmr_wme_gbmi$Estimate[1] * log(2) + qnorm(0.975) * mvmr_wme_gbmi$StdError[1] * log(2))

results_estimates$ltl_beta[4] <- mvmr_wmo_gbmi$Estimate[2]
results_estimates$ltl_p_value[4] <- mvmr_wmo_gbmi$Pvalue[2]
results_estimates$ltl_odds_ratio[4] <- exp(mvmr_wmo_gbmi$Estimate[2])
results_estimates$ltl_or_lower[4] <- exp(mvmr_wmo_gbmi$CILower[2])
results_estimates$ltl_or_upper[4] <- exp(mvmr_wmo_gbmi$CIUpper[2])

results_estimates$chip_beta[4] <- mvmr_wmo_gbmi$Estimate[1]
results_estimates$chip_p_value[4] <- mvmr_wmo_gbmi$Pvalue[1]
results_estimates$chip_odds_ratio[4] <- exp(mvmr_wmo_gbmi$Estimate[1] * log(2))
results_estimates$chip_or_lower[4] <- exp(mvmr_wmo_gbmi$CILower[1] * log(2))
results_estimates$chip_or_upper[4] <- exp(mvmr_wmo_gbmi$CIUpper[1] * log(2))

results_pleiotropy$heter_stat[1] <- mvmr_ivw_gbmi$Heter.Stat[1]
results_pleiotropy$heter_p[1] <- mvmr_ivw_gbmi$Heter.Stat[2]
results_pleiotropy$mr_presso_rss[1] <- mr_presso_res_gbmi$`MR-PRESSO results`$`Global Test`$RSSobs
results_pleiotropy$mr_presso_global_test_p[1] <- mr_presso_res_gbmi$`MR-PRESSO results`$`Global Test`$Pvalue
if (is.na(mr_presso_res_gbmi$`Main MR results`$`Causal Estimate`[3])) {
    results_pleiotropy$mr_presso_num_outliers[1] <- 0
}
results_pleiotropy$egger_intercept[1] <- mvmr_egger_gbmi$Intercept
results_pleiotropy$egger_intercept_p[1] <- mvmr_egger_gbmi$Pvalue.Int

results_instruments$total_snps[1] <- length(gbmi_out_dat$SNP)
results_instruments$ltl_snps[1] <- length(ltl_instruments_filtered_gbmi$SNP)
results_instruments$chip_snps[1] <- length(chip_instruments_filtered_gbmi$SNP)
results_instruments$overlap[1] <- length(instruments_intersection_gbmi)

results_instrument_strength$ltl_cond_f[1] <- as.numeric(sre_gbmi$exposure2[1])
results_instrument_strength$chip_cond_f[1] <- as.numeric(sre_gbmi$exposure1[1])

tsmr_res_gbmi <- mv_multiple(mvdat_gbmi, plots = T)

chip_plot_gbmi <- tsmr_res_gbmi$plot[[1]] +
    theme_bw() +
    ylab("Marginal SNP effect on IPF (GBMI meta-analysis)") +
    scale_x_continuous(
        limits = c(0, 0.2),
        breaks = seq(0, 0.2, by = 0.1)
    ) +
    scale_y_continuous(
        limits = c(-0.3, 0.1),
        breaks = c(-0.3, -0.2, -0.1, 0, 0.1)
    ) +
    theme(
        text = element_text(family = "alegreya_sans", size = 90),
    )
chip_plot_gbmi$layers <- append(
    chip_plot_gbmi$layers,
    list(
        geom_hline(yintercept = 0, colour = "darkgray"),
        geom_vline(xintercept = 0, colour = "darkgray")
    ),
    after = 0
)

ltl_plot_gbmi <- tsmr_res_gbmi$plot[[2]] +
    theme_bw() +
    ylab("Marginal SNP effect on IPF (GBMI meta-analysis)") +
    scale_x_continuous(
        limits = c(0, 0.2),
        breaks = seq(0, 0.2, by = 0.1)
    ) +
    scale_y_continuous(
        limits = c(-0.3, 0.1),
        breaks = c(-0.3, -0.2, -0.1, 0, 0.1)
    ) +
    theme(
        text = element_text(family = "alegreya_sans", size = 90),
    )
ltl_plot_gbmi$layers <- append(
    ltl_plot_gbmi$layers,
    list(
        geom_hline(yintercept = 0, colour = "darkgray"),
        geom_vline(xintercept = 0, colour = "darkgray")
    ),
    after = 0
)

chip_plot_gbmi <- chip_plot_gbmi + labs(y = NULL)
ltl_plot_gbmi <- ltl_plot_gbmi + labs(y = NULL)
gbmi_plot <- ltl_plot_gbmi + chip_plot_gbmi

gbmi_plot <- wrap_elements(gbmi_plot) +
    labs(tag = "Marginal SNP effect on IPF (GBMI meta-analysis)") +
    theme(
        plot.tag = element_text(family = "alegreya_sans", size = 90, angle = 90),
        plot.tag.position = "left"
    )

mr_presso_res_fg <- mr_presso(BetaOutcome = "betaYG", BetaExposure = c("betaX1", "betaX2"),
                              SdOutcome = "sebetaYG", SdExposure = c("sebetaX1", "sebetaX2"),
                              OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = F_data_fg,
                              NbDistribution = 1000, SignifThreshold = 0.05)

mvmr_pleio_res_fg <- pleiotropy_mvmr(r_input = F_data_fg, gencov = 0)
sre_fg <- strength_mvmr(r_input = F_data_fg, gencov = 0)

mvmr_ivw_fg <- mr_mvivw(mvinput_fg)
mvmr_egger_fg <- mr_mvegger(mvinput_fg)
mvmr_wme_fg <- mr_mvmedian(mvinput_fg)
mvmr_wmo_fg <- mv_mrmode(Bout = mvdat_fg$outcome_beta, 
                         Bexp = mvdat_fg$exposure_beta,
                         SEout = mvdat_fg$outcome_se, 
                         SEexp = mvdat_fg$exposure_se, 
                         Psi=0, CIMin = NA, CIMax = NA, CIStep = 0.001, alpha = 0.05,
                         residual="IVW", Mode="CM", weighting = "weighted", stderror = "simple", 
                         phi = 1, distribution = "normal", iterations = 10000)

ltl_instruments_filtered_fg <- ltl_instruments[ltl_instruments$SNP %in% fg_out_dat$SNP, ]
chip_instruments_filtered_fg <- chip_instruments[chip_instruments$SNP %in% fg_out_dat$SNP, ]
instruments_intersection_fg <- intersect(ltl_instruments_filtered_fg$SNP, chip_instruments_filtered_fg$SNP)

results_estimates$ltl_beta[5] <- mvmr_ivw_fg$Estimate[2]
results_estimates$ltl_se[5] <- mvmr_ivw_fg$StdError[2]
results_estimates$ltl_p_value[5] <- mvmr_ivw_fg$Pvalue[2]
results_estimates$ltl_odds_ratio[5] <- exp(mvmr_ivw_fg$Estimate[2])
results_estimates$ltl_or_lower[5] <- exp(mvmr_ivw_fg$Estimate[2] - qnorm(0.975) * mvmr_ivw_fg$StdError[2])
results_estimates$ltl_or_upper[5] <- exp(mvmr_ivw_fg$Estimate[2] + qnorm(0.975) * mvmr_ivw_fg$StdError[2])

results_estimates$chip_beta[5] <- mvmr_ivw_fg$Estimate[1]
results_estimates$chip_se[5] <- mvmr_ivw_fg$StdError[1]
results_estimates$chip_p_value[5] <- mvmr_ivw_fg$Pvalue[1]
results_estimates$chip_odds_ratio[5] <- exp(mvmr_ivw_fg$Estimate[1] * log(2))
results_estimates$chip_or_lower[5] <- exp(mvmr_ivw_fg$Estimate[1] * log(2) - qnorm(0.975) * mvmr_ivw_fg$StdError[1] * log(2))
results_estimates$chip_or_upper[5] <- exp(mvmr_ivw_fg$Estimate[1] * log(2) + qnorm(0.975) * mvmr_ivw_fg$StdError[1] * log(2))

results_estimates$ltl_beta[6] <- mvmr_egger_fg$Estimate[2]
results_estimates$ltl_se[6] <- mvmr_egger_fg$StdError.Est[2]
results_estimates$ltl_p_value[6] <- mvmr_egger_fg$Pvalue.Est[2]
results_estimates$ltl_odds_ratio[6] <- exp(mvmr_egger_fg$Estimate[2])
results_estimates$ltl_or_lower[6] <- exp(mvmr_egger_fg$Estimate[2] - qnorm(0.975) * mvmr_egger_fg$StdError.Est[2])
results_estimates$ltl_or_upper[6] <- exp(mvmr_egger_fg$Estimate[2] + qnorm(0.975) * mvmr_egger_fg$StdError.Est[2])

results_estimates$chip_beta[6] <- mvmr_egger_fg$Estimate[1]
results_estimates$chip_se[6] <- mvmr_egger_fg$StdError.Est[1]
results_estimates$chip_p_value[6] <- mvmr_egger_fg$Pvalue.Est[1]
results_estimates$chip_odds_ratio[6] <- exp(mvmr_egger_fg$Estimate[1] * log(2))
results_estimates$chip_or_lower[6] <- exp(mvmr_egger_fg$Estimate[1] * log(2) - qnorm(0.975) * mvmr_egger_fg$StdError.Est[1] * log(2))
results_estimates$chip_or_upper[6] <- exp(mvmr_egger_fg$Estimate[1] * log(2) + qnorm(0.975) * mvmr_egger_fg$StdError.Est[1] * log(2))

results_estimates$ltl_beta[7] <- mvmr_wme_fg$Estimate[2]
results_estimates$ltl_se[7] <- mvmr_wme_fg$StdError[2]
results_estimates$ltl_p_value[7] <- mvmr_wme_fg$Pvalue[2]
results_estimates$ltl_odds_ratio[7] <- exp(mvmr_wme_fg$Estimate[2])
results_estimates$ltl_or_lower[7] <- exp(mvmr_wme_fg$Estimate[2] - qnorm(0.975) * mvmr_ivw_fg$StdError[2])
results_estimates$ltl_or_upper[7] <- exp(mvmr_wme_fg$Estimate[2] + qnorm(0.975) * mvmr_ivw_fg$StdError[2])

results_estimates$chip_beta[7] <- mvmr_wme_fg$Estimate[1]
results_estimates$chip_se[7] <- mvmr_wme_fg$StdError[1]
results_estimates$chip_p_value[7] <- mvmr_wme_fg$Pvalue[1]
results_estimates$chip_odds_ratio[7] <- exp(mvmr_wme_fg$Estimate[1] * log(2))
results_estimates$chip_or_lower[7] <- exp(mvmr_wme_fg$Estimate[1] * log(2) - qnorm(0.975) * mvmr_wme_fg$StdError[1] * log(2))
results_estimates$chip_or_upper[7] <- exp(mvmr_wme_fg$Estimate[1] * log(2) + qnorm(0.975) * mvmr_wme_fg$StdError[1] * log(2))

results_estimates$ltl_beta[8] <- mvmr_wmo_fg$Estimate[2]
results_estimates$ltl_p_value[8] <- mvmr_wmo_fg$Pvalue[2]
results_estimates$ltl_odds_ratio[8] <- exp(mvmr_wmo_fg$Estimate[2])
results_estimates$ltl_or_lower[8] <- exp(mvmr_wmo_fg$CILower[2])
results_estimates$ltl_or_upper[8] <- exp(mvmr_wmo_fg$CIUpper[2])

results_estimates$chip_beta[8] <- mvmr_wmo_fg$Estimate[1]
results_estimates$chip_p_value[8] <- mvmr_wmo_fg$Pvalue[1]
results_estimates$chip_odds_ratio[8] <- exp(mvmr_wmo_fg$Estimate[1] * log(2))
results_estimates$chip_or_lower[8] <- exp(mvmr_wmo_fg$CILower[1] * log(2))
results_estimates$chip_or_upper[8] <- exp(mvmr_wmo_fg$CIUpper[1] * log(2))

results_pleiotropy$heter_stat[2] <- mvmr_ivw_fg$Heter.Stat[1]
results_pleiotropy$heter_p[2] <- mvmr_ivw_fg$Heter.Stat[2]
results_pleiotropy$mr_presso_rss[2] <- mr_presso_res_fg$`MR-PRESSO results`$`Global Test`$RSSobs
results_pleiotropy$mr_presso_global_test_p[2] <- mr_presso_res_fg$`MR-PRESSO results`$`Global Test`$Pvalue
if (is.na(mr_presso_res_fg$`Main MR results`$`Causal Estimate`[3])) {
    results_pleiotropy$mr_presso_num_outliers[2] <- 0
}
results_pleiotropy$egger_intercept[2] <- mvmr_egger_fg$Intercept
results_pleiotropy$egger_intercept_p[2] <- mvmr_egger_fg$Pvalue.Int

results_instruments$total_snps[2] <- length(fg_out_dat$SNP)
results_instruments$ltl_snps[2] <- length(ltl_instruments_filtered_fg$SNP)
results_instruments$chip_snps[2] <- length(chip_instruments_filtered_fg$SNP)
results_instruments$overlap[2] <- length(instruments_intersection_fg)

results_instrument_strength$ltl_cond_f[2] <- as.numeric(sre_fg$exposure2[1])
results_instrument_strength$chip_cond_f[2] <- as.numeric(sre_fg$exposure1[1])

tsmr_res_fg <- mv_multiple(mvdat_fg, plots = T)

chip_plot_fg <- tsmr_res_fg$plot[[1]] +
    theme_bw() +
    ylab("Marginal SNP effect on IPF (FinnGen GWAS)") +
    scale_x_continuous(
        limits = c(0, 0.2),
        breaks = seq(0, 0.2, by = 0.1)
    ) +
    scale_y_continuous(
        limits = c(-0.3, 0.1),
        breaks = c(-0.3, -0.2, -0.1, 0, 0.1)
    ) +
    theme(
        text = element_text(family = "alegreya_sans", size = 90),
    )
chip_plot_fg$layers <- append(
    chip_plot_fg$layers,
    list(
        geom_hline(yintercept = 0, colour = "darkgray"),
        geom_vline(xintercept = 0, colour = "darkgray")
    ),
    after = 0
)

ltl_plot_fg <- tsmr_res_fg$plot[[2]] +
    theme_bw() +
    ylab("Marginal SNP effect on IPF (FinnGen GWAS)") +
    scale_x_continuous(
        limits = c(0, 0.2),
        breaks = seq(0, 0.2, by = 0.1)
    ) +
    scale_y_continuous(
        limits = c(-0.3, 0.1),
        breaks = c(-0.3, -0.2, -0.1, 0, 0.1)
    ) +
    theme(
        text = element_text(family = "alegreya_sans", size = 90),
    )
ltl_plot_fg$layers <- append(
    ltl_plot_fg$layers,
    list(
        geom_hline(yintercept = 0, colour = "darkgray"),
        geom_vline(xintercept = 0, colour = "darkgray")
    ),
    after = 0
)

chip_plot_fg <- chip_plot_fg + labs(y = NULL)
ltl_plot_fg <- ltl_plot_fg + labs(y = NULL)
fg_plot <- ltl_plot_fg + chip_plot_fg

fg_plot <- wrap_elements(fg_plot) +
    labs(tag = "Marginal SNP effect on IPF (FinnGen GWAS)") +
    theme(
        plot.tag = element_text(family = "alegreya_sans", size = 90, angle = 90),
        plot.tag.position = "left"
    )

write.table(results_estimates, "results/ltl_chip_ipf_mvmr_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(results_pleiotropy, "results/ltl_chip_ipf_mvmr_results_pleiotropy.txt", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(results_instruments, "results/ltl_chip_ipf_mvmr_results_instruments.txt", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(results_instrument_strength, "results/ltl_chip_ipf_mvmr_results_instrument_strength.txt", quote = FALSE, row.names = FALSE, sep = "\t")

final_plot <- gbmi_plot / fg_plot

ggsave("plots/mvmr_plot.png", final_plot, width = 9, height = 12, dpi = 600)
