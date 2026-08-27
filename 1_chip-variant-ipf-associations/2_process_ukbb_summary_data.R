raw_ipf_data <- read.table(
    "data/filtered_ipf_summary_stats.txt",
    header = TRUE)

raw_ltl_data <- read.table(
    "data/filtered_ltl_summary_stats.txt",
    header = TRUE)

raw_published_data <- read.csv(
    "data/published_ch_candidate_mutations.csv",
    header = TRUE)

ids_df <- data.frame(
    ID = raw_published_data$ID
)

filtered_ipf_data <- merge(raw_ipf_data, ids_df, by = "ID")
filtered_ipf_data <- filtered_ipf_data[!duplicated(filtered_ipf_data), ]

filtered_ltl_data <- merge(raw_ltl_data, ids_df, by = "ID")
filtered_ltl_data <- filtered_ltl_data[!duplicated(filtered_ltl_data), ]

rm(ids_df, raw_ipf_data, raw_ltl_data)

ipf_data_to_join <- data.frame(
    ID = filtered_ipf_data$ID,
    IPF_N = filtered_ipf_data$N,
    IPF_freq = filtered_ipf_data$A1FREQ,
    IPF_alt_count = round(2 * filtered_ipf_data$N * filtered_ipf_data$A1FREQ),
    IPF_beta = filtered_ipf_data$BETA,
    IPF_odds_ratio = exp(filtered_ipf_data$BETA),
    IPF_se = filtered_ipf_data$SE,
    IPF_log10p = filtered_ipf_data$LOG10P,
    IPF_p_value = 10 ^ (-1 * filtered_ipf_data$LOG10P)
)

ltl_data_to_join <- data.frame(
    ID = filtered_ltl_data$ID,
    LTL_N = filtered_ltl_data$N,
    LTL_freq = filtered_ltl_data$A1FREQ,
    LTL_alt_count = round(2 * filtered_ltl_data$N * filtered_ltl_data$A1FREQ),
    LTL_beta = filtered_ltl_data$BETA,
    LTL_se = filtered_ltl_data$SE,
    LTL_log10p = filtered_ltl_data$LOG10P,
    LTL_p_value = 10 ^ (-1 * filtered_ltl_data$LOG10P)
)

rm(filtered_ipf_data, filtered_ltl_data)

ukbb_data_to_join <- merge(ipf_data_to_join, ltl_data_to_join, by = "ID",
                           all = TRUE)

rm(ipf_data_to_join, ltl_data_to_join)

published_coords <- strsplit(raw_published_data$ID, ":", fixed = TRUE)
published_chromosomes <- paste("chr", sapply(published_coords,
                                             function(x) x[1]),
                               sep = "")
published_positions <-  as.integer(sapply(published_coords, function(x) x[2]))
published_refs <- sapply(published_coords, function(x) x[3])
published_alts <- sapply(published_coords, function(x) x[4])

rm(published_coords)

published_data_to_join <- data.frame(
    GENE = raw_published_data$GENE,
    ID = raw_published_data$ID,
    CHROM = published_chromosomes,
    POS = published_positions,
    REF = published_refs,
    ALT = published_alts,
    BickMeanVAF = raw_published_data$BickMeanVAF,
    GutierrezRodriguesMeanVAF = raw_published_data$GutierrezRodriguesMeanVAF,
    KarMeanVAF = raw_published_data$KarMeanVAF
)

rm(published_chromosomes, published_positions, published_refs, published_alts,
   raw_published_data)

final_table <- merge(published_data_to_join, ukbb_data_to_join, by = "ID")

write.csv(final_table, "results/ukbb_association_statistics_for_ch_mutations.csv",
          quote = FALSE, row.names = FALSE)
