# ============================================================================
# Build proliferation-related combined GMT file (curated 18 gene sets)
# Sources: Hallmark + KEGG + PROGENy
# Only core proliferation & growth signaling pathways
# ============================================================================

library(msigdbr)
library(dplyr)
library(scMMR)

gmt_dir <- system.file("extdata", "gmt", package = "scMMR")

write_gmt <- function(gene_sets, file) {
  con <- file(file, "w")
  on.exit(close(con))
  for (name in names(gene_sets)) {
    genes <- unique(gene_sets[[name]])
    line <- paste(c(name, "", genes), collapse = "\t")
    writeLines(line, con)
  }
}


# ── 1. Hallmark proliferation pathways (9 sets) ─────────────────────────────
cat("=== 1. Hallmark proliferation pathways ===\n")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol)

prolif_hallmark_names <- c(
  # Core cell cycle
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  # Growth signaling
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

prolif_hallmark <- hallmark %>% filter(gs_name %in% prolif_hallmark_names)
hallmark_list <- split(prolif_hallmark[["gene_symbol"]], prolif_hallmark[["gs_name"]])
hallmark_list <- lapply(hallmark_list, unique)

cat("Hallmark sets:", length(hallmark_list), "\n")
for (n in names(hallmark_list)) {
  cat("  ", n, ":", length(hallmark_list[[n]]), "genes\n")
}


# ── 2. KEGG proliferation pathways (5 sets) ──────────────────────────────────
cat("\n=== 2. KEGG proliferation pathways ===\n")

kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  select(gs_name, gene_symbol)

prolif_kegg_names <- c(
  "KEGG_CELL_CYCLE",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_ERBB_SIGNALING_PATHWAY",
  "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_WNT_SIGNALING_PATHWAY"
)

prolif_kegg <- kegg %>% filter(gs_name %in% prolif_kegg_names)
kegg_list <- split(prolif_kegg[["gene_symbol"]], prolif_kegg[["gs_name"]])
kegg_list <- lapply(kegg_list, unique)

cat("KEGG sets:", length(kegg_list), "\n")
for (n in names(kegg_list)) {
  cat("  ", n, ":", length(kegg_list[[n]]), "genes\n")
}


# ── 3. PROGENy proliferation pathways (4 sets) ──────────────────────────────
cat("\n=== 3. PROGENy proliferation pathways ===\n")

progeny_file <- file.path(gmt_dir, "progeny.human.top500.gmt")
progeny_all <- read_gmt(progeny_file)

progeny_prolif_names <- c("EGFR", "MAPK", "PI3K", "WNT")
progeny_list <- progeny_all[names(progeny_all) %in% progeny_prolif_names]
names(progeny_list) <- paste0("PROGENY_", names(progeny_list))

cat("PROGENy sets:", length(progeny_list), "\n")
for (n in names(progeny_list)) {
  cat("  ", n, ":", length(progeny_list[[n]]), "genes\n")
}


# ── 4. Combine and save ─────────────────────────────────────────────────────
cat("\n=== 4. Combining ===\n")

all_prolif <- c(hallmark_list, kegg_list, progeny_list)
cat("Total gene sets:", length(all_prolif), "\n")
cat("Total unique genes:", length(unique(unlist(all_prolif))), "\n")

out_file <- "/home/oyh/project/para/python_pacakge/scMMR/inst/extdata/gmt/proliferation.combined.gmt"
write_gmt(all_prolif, out_file)
cat("Written:", out_file, "\n")


# ── 5. Summary ──────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat(sprintf("%-45s %5s  %s\n", "Pathway", "Genes", "Source"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (n in names(all_prolif)) {
  src <- ifelse(grepl("^HALLMARK_", n), "Hallmark",
         ifelse(grepl("^KEGG_", n), "KEGG",
         ifelse(grepl("^PROGENY_", n), "PROGENy", "Other")))
  cat(sprintf("%-45s %5d  %s\n", n, length(all_prolif[[n]]), src))
}
