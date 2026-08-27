library(easylift)
library(enshuman)
library(GenomicRanges)
library(topr)

tert_transcript <- read.csv(
    "data/ensembl_tert_transcript.csv",
    header = TRUE) # Downloadable from: 
                   # https://jun2026.archive.ensembl.org/Homo_sapiens/Transcript/
                   # Exons?db=core;g=ENSG00000164362;r=5:1253147-1295086;t=ENST00000310581
                   # (click on the Excel symbol in the top right had corner of the table
                   # and then "Download what you see")

tert_start <- as.integer(
    gsub(",", "",
         tert_transcript$Start[!is.na(tert_transcript$No.)
                               & tert_transcript$No. == 1]))

tert_exon_1 <- tert_transcript$Sequence[!is.na(tert_transcript$No.)
                                        & tert_transcript$No. == 1]

tert_exon_1 <- gsub("[[:space:]]", "", tert_exon_1)

atg_pos <- as.integer(regexpr("ATG", tert_exon_1))

tert_start_codon <- tert_start - atg_pos + 1

# According to https://doi.org/10.3389/fcell.2023.1286683:
# "The TERT gene is located approximately one megabase from the
# end of the short arm of chromosome 5, and its core promoter is
# delimited by nucleotides located at positions −330 upstream and
# +37 downstream of the TERT ATG start site"

tert_promoter_bounds <- tert_start_codon + c(-37, 330)

rm(atg_pos, tert_exon_1, tert_start, tert_start_codon, tert_transcript)

gr_st1 <- read.csv(
    paste("data/gutierrez-rodrigues_2024_ST1-",
          "Clinical-characteristics-and-genomic-data-of-207-TBD-patients.csv",
          sep = ""),
    header = TRUE)

gr_st1 <- gr_st1[gr_st1$Genomic.Location != "",]

gr_locations <- strsplit(gr_st1$Genomic.Location, ":", fixed = TRUE)

gr_coords <- data.frame(
  chr = sapply(gr_locations, function(x) x[1]),
  pos = as.integer(sapply(gr_locations, function(x) x[2]))
)

download.file(
    paste("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
          "hg19ToHg38.over.chain.gz", sep = ""),
    "hg19ToHg38.over.chain.gz"
)

gr_hg19_coordinates <- GRanges(
    seqnames = gr_coords$chr,
    ranges = IRanges(start = gr_coords$pos, width = 1)
)
genome(gr_hg19_coordinates) <- "hg19"

gr_hg38_coordinates <- easylift(gr_hg19_coordinates, "hg38",
                                "hg19ToHg38.over.chain.gz")

gr_chromosomes <- sub("chr", "", as.character(seqnames(gr_hg38_coordinates)))
gr_positions <- start(gr_hg38_coordinates)

file.remove("hg19ToHg38.over.chain.gz")

gr_bases <- strsplit(gr_st1$Ref.Alt.Allele, " / ", fixed = TRUE)
gr_ref_alleles <- sapply(gr_bases, function(x) x[1])
gr_alt_alleles <- sapply(gr_bases, function(x) x[2])

gr_identifiers <- paste(gr_chromosomes, gr_positions,
                        gr_ref_alleles, gr_alt_alleles, sep = ":")

gr_somatic_variants <- data.frame(
    pt = gr_st1$Study.ID[gr_st1$Classification == "Somatic"],
    vaf = as.numeric(gr_st1$Allele.Frequency[
        gr_st1$Classification == "Somatic"]),
    id = gr_identifiers[gr_st1$Classification == "Somatic"]
)

gr_somatic_variants$mean_vaf <- ave(gr_somatic_variants$vaf,
                                    gr_somatic_variants$id,
                                    FUN = mean)

gr_variants <- data.frame(
    ID = gr_somatic_variants$id,
    GutierrezRodriguesMeanVAF = gr_somatic_variants$mean_vaf
)

gr_variants <- gr_variants[!duplicated(gr_variants), ]

gr_temp_vars <- ls(pattern = "^gr_")
gr_temp_vars <- gr_temp_vars[gr_temp_vars != "gr_variants"]
rm(list = gr_temp_vars)
rm(gr_temp_vars)

bick_st3 <- read.csv(
    paste("data/bick_2020_ST3-",
          "CHIP-mutations-identified-in-TOPMed-participants.csv",
          sep = ""),
    header = TRUE)

bick_variants <- data.frame(
    ID = paste(bick_st3$CHROM, bick_st3$POS, bick_st3$REF, bick_st3$ALT,
               sep = ":"),
    VAF = bick_st3$VAF
)

bick_variants$BickMeanVAF <- ave(bick_variants$VAF,
                                 bick_variants$ID,
                                 FUN = mean)

bick_variants$VAF <- NULL

bick_variants <- bick_variants[!duplicated(bick_variants), ]

rm(bick_st3)

kar_st4 <- read.csv(
    paste("data/kar_2022_ST4-",
          "List-of-somatic-driver-mutations-identified-in-the-study-cohort",
          ".csv",
          sep = ""),
    header = TRUE)

kar_variants <- data.frame(
    ID = paste(sub("^chr", "", kar_st4$CHROM),
               kar_st4$POS, kar_st4$REF, kar_st4$ALT,
               sep = ":"),
    VAF = kar_st4$VAF
)

kar_variants$KarMeanVAF <- ave(kar_variants$VAF,
                               kar_variants$ID,
                               FUN = mean)

kar_variants$VAF <- NULL

kar_variants <- kar_variants[!duplicated(kar_variants), ]

rm(kar_st4)

published_ch_candidate_mutations <- bick_variants |>
    merge(gr_variants, by = "ID", all = TRUE) |>
    merge(kar_variants, by = "ID", all = TRUE)

rm(bick_variants, gr_variants, kar_variants)

published_ch_candidate_mutations <- cbind(GENE = NA,
                                          published_ch_candidate_mutations)

for (i in 1:nrow(published_ch_candidate_mutations)) {
    i_coords <- strsplit(published_ch_candidate_mutations$ID[i], ":")[[1]]
    i_chrom <- as.character(i_coords[1])
    i_pos <- as.integer(i_coords[2])
    if (i_chrom == 5
        && i_pos >= tert_promoter_bounds[1]
        && i_pos <= tert_promoter_bounds[2]) {
        published_ch_candidate_mutations$GENE[i] <- "TERTp"
    } else {
        i_gene_list <- hg38[hg38$chrom == i_chrom
                            & hg38$gene_start <= i_pos
                            & hg38$gene_end >= i_pos, ]
        if (length(i_gene_list$gene_symbol) > 1) {
            published_ch_candidate_mutations$GENE[
                i] <- annotate_with_nearest_gene(data.frame(
                    CHROM = paste("chr", i_chrom, sep = ""),
                    POS = i_pos
                ))$Gene_Symbol
        } else {
            published_ch_candidate_mutations$GENE[i] <- i_gene_list$gene_symbol
        }
    }
}

rm(i, i_chrom, i_coords, i_gene_list, i_pos)

write.csv(published_ch_candidate_mutations,
          "data/published_ch_candidate_mutations.csv",
          quote = FALSE, row.names = FALSE)

gene_list <- data.frame(
    GENE = unique(published_ch_candidate_mutations$GENE)
)

gene_list$CHROM <- NA
gene_list$START <- NA
gene_list$END <- NA

for (i in 1:nrow(gene_list)) {
    if (gene_list$GENE[i] == "TERTp") {
        gene_list$CHROM[i] <- 5
        gene_list$START[i] <- tert_promoter_bounds[1]
        gene_list$END[i] <- tert_promoter_bounds[2]
    } else {
        i_gene_coords <- get_gene_coords(gene_list$GENE[i])
        gene_list$CHROM[i] <- i_gene_coords$chrom
        gene_list$START[i] <- i_gene_coords$gene_start
        gene_list$END[i] <- i_gene_coords$gene_end
    }
}

write.table(gene_list, "data/published_ch_candiate_mutations_gene_list.txt", quote = FALSE, row.names = FALSE)
