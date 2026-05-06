# =============================================================================
# Multi-omics analysis reveals a non-coding variant
# connecting CRHR1 signalling to coronary artery disease
#
# GWAS variant prioritisation via sequential epigenomic integration:
#   1. hg18 -> hg19 liftOver
#   2. p < 0.05 filtering
#   3. DNase overlap
#   4. H3K27ac overlap
#   5. H3K36me3 overlap
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

# -----------------------------------------------------------------------------
# 1. GWAS liftOver: hg18 -> hg19
# -----------------------------------------------------------------------------

gwas <- fread("GWAS.csv")

chr_pos <- tstrsplit(gwas$`chr_pos_(b36)`, ":", fixed = TRUE)
gwas[, chr_hg18 := chr_pos[[1]]]
gwas[, pos_hg18 := as.integer(chr_pos[[2]])]

gwas <- gwas[!is.na(chr_hg18) & !is.na(pos_hg18)]

gr_hg18 <- GRanges(
  seqnames = gwas$chr_hg18,
  ranges   = IRanges(start = gwas$pos_hg18, end = gwas$pos_hg18),
  SNP      = gwas$SNP
)

chain <- import.chain("hg18ToHg19.over.chain.gz")
gr_hg19 <- unlist(liftOver(gr_hg18, chain))

lifted <- data.table(
  SNP      = mcols(gr_hg19)$SNP,
  chr_hg19 = as.character(seqnames(gr_hg19)),
  pos_hg19 = start(gr_hg19)
)

gwas_hg19 <- merge(gwas, lifted, by = "SNP", all.x = TRUE, sort = FALSE)

mapped   <- gwas_hg19[!is.na(chr_hg19) & !is.na(pos_hg19)]
unmapped <- gwas_hg19[ is.na(chr_hg19) | is.na(pos_hg19)]

fwrite(mapped,   "gwas_hg19_mapped.csv")
fwrite(unmapped, "gwas_hg19_unmapped.csv")

cat(sprintf("%-24s: %d\n", "Total GWAS variants", nrow(gwas_hg19)))
cat(sprintf("%-24s: %d (%.2f%%)\n",
            "Successfully mapped",
            nrow(mapped),
            100 * nrow(mapped) / nrow(gwas_hg19)))


# -----------------------------------------------------------------------------
# 2. Sequential epigenomic filtering
# -----------------------------------------------------------------------------

# Step 1 — p < 0.05
gwas_p005 <- mapped[!is.na(pvalue) & pvalue < 0.05]
cat(sprintf("%-24s: %d\n", "After p < 0.05", nrow(gwas_p005)))

gr_gwas <- GRanges(
  seqnames = gwas_p005$chr_hg19,
  ranges   = IRanges(gwas_p005$pos_hg19, gwas_p005$pos_hg19),
  SNP      = gwas_p005$SNP
)


# Step 2 — DNase-seq
dnase <- fread("dna_seq.bed", header = FALSE, sep = "\t", fill = TRUE)
dnase <- dnase[, .(chr = V1, start = as.integer(V2), end = as.integer(V3))]
dnase <- dnase[!is.na(start) & !is.na(end) & end > start]

gr_dnase <- GRanges(
  seqnames = dnase$chr,
  ranges   = IRanges(start = dnase$start + 1L, end = dnase$end)
)

hits_dnase <- findOverlaps(gr_gwas, gr_dnase)
gwas_dnase <- gwas_p005[unique(queryHits(hits_dnase))]
cat(sprintf("%-24s: %d\n", "After DNase", nrow(gwas_dnase)))


# Step 3 — H3K27ac (merged peaks)
h3k27ac <- fread("GSM906392_UCSD.Aorta.H3K27ac.STL003.bed",
                 header = FALSE, sep = "\t")
setnames(h3k27ac, c("chr","start","end","name","score","strand"))

gr_h3k27 <- reduce(GRanges(
  seqnames = h3k27ac$chr,
  ranges   = IRanges(start = h3k27ac$start + 1L, end = h3k27ac$end)
))

gr_after_dnase <- GRanges(
  seqnames = gwas_dnase$chr_hg19,
  ranges   = IRanges(gwas_dnase$pos_hg19, gwas_dnase$pos_hg19),
  SNP      = gwas_dnase$SNP
)

hits_h3k27 <- findOverlaps(gr_after_dnase, gr_h3k27)
gwas_h3k27 <- gwas_dnase[unique(queryHits(hits_h3k27))]
cat(sprintf("%-24s: %d\n", "After H3K27ac", nrow(gwas_h3k27)))


# Step 4 — H3K36me3 (merged peaks)
h3k36 <- fread("GSM910566_UCSD.Aorta.H3K36me3.STL003.bed",
               header = FALSE, sep = "\t")
setnames(h3k36, c("chr","start","end","name","score","strand"))

gr_h3k36 <- reduce(GRanges(
  seqnames = h3k36$chr,
  ranges   = IRanges(start = h3k36$start + 1L, end = h3k36$end)
))

gr_after_h3k27 <- GRanges(
  seqnames = gwas_h3k27$chr_hg19,
  ranges   = IRanges(gwas_h3k27$pos_hg19, gwas_h3k27$pos_hg19),
  SNP      = gwas_h3k27$SNP
)

hits_h3k36 <- findOverlaps(gr_after_h3k27, gr_h3k36)
final_set  <- gwas_h3k27[unique(queryHits(hits_h3k36))]
final_set  <- unique(final_set, by = "SNP")

cat(sprintf("%-24s: %d\n", "After H3K36me3", nrow(final_set)))


# -----------------------------------------------------------------------------
# 3. Output
# -----------------------------------------------------------------------------

fwrite(final_set, "CAD_p005_DNase_H3K27ac_H3K36me3.csv")

cat("\nOutput written to:\n")
cat("  CAD_p005_DNase_H3K27ac_H3K36me3.csv\n")