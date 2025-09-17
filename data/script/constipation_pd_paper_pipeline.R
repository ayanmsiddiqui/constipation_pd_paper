# constipation_pd_paper_pipeline.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tibble)
  library(data.table)
  library(vegan)
  library(phyloseq)
  library(limma)
  library(edgeR)
  library(patchwork)
  library(scales)
  library(stringr)
  library(mclust)
  library(cluster)
  library(ggrepel)
  library(forcats)
})

# ---------- misc helpers ----------
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse=" ")))
safe_mkdir <- function(path) if (!dir.exists(path)) dir.create(path, recursive=TRUE, showWarnings=FALSE)
fmt_p <- function(x) {
  y <- suppressWarnings(as.numeric(x)); y[is.na(y)] <- NA_real_; y <- pmin(pmax(y,0),1)
  ifelse(y < 1e-16, "<1e-16", ifelse(y < 1e-3, formatC(y, format="e", digits=2), signif(y, 3)))
}
save_figure <- function(plot_obj, file_base, width=7.2, height=4.5, dpi=400) {
  ggsave(paste0(file_base, ".png"), plot_obj, width=width, height=height, dpi=dpi, bg="white")
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(paste0(file_base, ".pdf"), plot_obj, width=width, height=height, device=dev)
}
# robust path finder: tries root_dir/path then root_dir/data/basename(path)
locate <- function(root_dir, path) {
  p1 <- file.path(root_dir, path); if (file.exists(p1)) return(p1)
  p2 <- file.path(root_dir, "data", basename(path)); if (file.exists(p2)) return(p2)
  stop(sprintf("File '%s' not found in '%s' or '%s/data'.", path, root_dir, root_dir))
}

# ---------- simple HTML export ----------
htmlEscape <- function(s) { s <- as.character(s); s <- gsub("&","&amp;",s,fixed=TRUE); s <- gsub("<","&lt;",s,fixed=TRUE); s <- gsub(">","&gt;",s,fixed=TRUE); s }
write_html_table <- function(df, file_html, title=NULL, footnotes=NULL) {
  css <- "<style>body{font-family:Arial}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left}th{background:#f5f5f5}</style>"
  con <- file(file_html, "w", encoding="UTF-8"); on.exit(close(con))
  writeLines(paste0("<html><head><meta charset='utf-8'>", css, "</head><body>"), con)
  if (!is.null(title)) writeLines(paste0("<h3>", htmlEscape(title), "</h3>"), con)
  writeLines("<table><thead><tr>", con); for (nm in names(df)) writeLines(paste0("<th>", htmlEscape(nm), "</th>"), con)
  writeLines("</tr></thead><tbody>", con)
  apply(df, 1, function(row){ writeLines("<tr>", con); for (v in row) writeLines(paste0("<td>", htmlEscape(v), "</td>"), con); writeLines("</tr>", con) })
  writeLines("</tbody>", con)
  if (!is.null(footnotes) && length(footnotes)) { writeLines("<tfoot><tr><td colspan='100%'><ul>", con); for (ft in footnotes) writeLines(paste0("<li>", htmlEscape(ft), "</li>"), con); writeLines("</ul></td></tr></tfoot>", con) }
  writeLines("</table></body></html>", con)
}
save_table_with_html <- function(df, csv_path, html_path, title=NULL, footnotes=character()) {
  fwrite(df, csv_path); write_html_table(df, html_path, title=title, footnotes=footnotes)
}

# ---------- plotting theme & scales ----------
base_font_size <- 11
theme_base <- function() {
  theme_classic(base_size=base_font_size) +
    theme(legend.position="top",
          legend.title=element_text(size=base_font_size-1),
          legend.text =element_text(size=base_font_size-1),
          plot.caption = element_text(size=base_font_size-2, hjust=0))
}
scale_fill_group <- function() scale_fill_manual(values=c(control="#4C78A8", case="#F58518"))
scale_colour_group <- function() scale_colour_manual(values=c(control="#4C78A8", case="#F58518"))

# ---------- labels / transforms ----------
short_tax_label <- function(x) {
  last <- vapply(strsplit(as.character(x), "\\|"), function(v) v[length(v)], "", USE.NAMES=FALSE)
  gsub("^[a-z]__", "", ifelse(nchar(last)>0, last, x))
}
short_path_label <- function(x, max_char=80) {
  s <- gsub("^PWY-\\d+[: ]*", "", as.character(x)); s <- gsub("_", " ", s); stringr::str_trunc(s, max_char)
}
normalize_rel <- function(ps) transform_sample_counts(ps, function(x) x/sum(x))
clr_transform <- function(mat, pseudo=1e-6) {
  m <- as.matrix(mat); m <- sweep(m, 2, colSums(m), "/"); m[!is.finite(m)] <- 0
  m <- m + pseudo
  log_m <- log(m); log_m - matrix(colMeans(log_m), nrow=nrow(log_m), ncol=ncol(log_m), byrow=TRUE)
}

# ---------- ID alignment ----------
.clean_ids <- function(x) { x <- as.character(x); x <- trimws(x); x <- gsub("\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$|_R[12]$|_R[12]_|-R[12]-","",x,ignore.case=TRUE); x <- gsub("[^A-Za-z0-9]+","_",x); tolower(x) }
align_counts_meta <- function(counts, meta_df, id_col) {
  stopifnot(id_col %in% names(meta_df))
  ids_counts_raw <- colnames(counts); ids_meta_raw <- meta_df[[id_col]]
  ids_counts_norm <- .clean_ids(ids_counts_raw); ids_meta_norm <- .clean_ids(ids_meta_raw)
  m <- match(ids_counts_norm, ids_meta_norm); has_match <- !is.na(m)
  counts2 <- counts[, has_match, drop=FALSE]; meta2 <- meta_df[m[has_match], , drop=FALSE]
  colnames(counts2) <- meta2[[id_col]]
  diag <- c(
    sprintf("Counts columns: %d; metadata rows: %d; matched: %d", length(ids_counts_raw), length(ids_meta_raw), sum(has_match)),
    if (any(!has_match)) paste("Unmatched count IDs:", paste(head(ids_counts_raw[!has_match], 10), collapse=", ")) else "Unmatched count IDs: none",
    {
      meta_only <- ids_meta_raw[setdiff(seq_along(ids_meta_raw), m[has_match])]
      if (length(meta_only)) paste("Unmatched meta IDs:", paste(head(meta_only, 10), collapse=", ")) else "Unmatched meta IDs: none"
    }
  )
  list(counts=counts2, meta=meta2, diagnostics=diag)
}

# ---------- group/covariates ----------
standardize_group_values <- function(vec, case_values) {
  v <- tolower(as.character(vec)); case_vals <- unique(tolower(c(case_values, "y","yes","true","1","t","constipated")))
  factor(ifelse(v %in% case_vals, "case", "control"), levels=c("control","case"))
}
prepare_design <- function(meta_df, group_col="Group", covars=character()) {
  md <- meta_df; stopifnot(group_col %in% names(md)); md[[group_col]] <- droplevels(factor(md[[group_col]]))
  keep_cov <- c()
  for (nm in covars) {
    if (!nm %in% names(md)) next
    v <- md[[nm]]
    if (is.character(v) || is.logical(v)) v <- factor(v)
    if (is.factor(v)) { v <- droplevels(v); if (nlevels(v) >= 2) keep_cov <- c(keep_cov, nm) }
    else if (is.numeric(v)) { if (sum(!is.na(v)) >= 2 && length(unique(v[!is.na(v)])) >= 2) keep_cov <- c(keep_cov, nm) }
    md[[nm]] <- v
  }
  list(meta=md, covars=keep_cov)
}
align_for_model <- function(mat, meta, id_col, covars) {
  sample_ids <- colnames(mat)
  meta2 <- meta %>% dplyr::filter(.data[[id_col]] %in% sample_ids) %>% dplyr::arrange(match(.data[[id_col]], sample_ids))
  design_vars <- c("Group", covars); cc <- stats::complete.cases(meta2[, design_vars, drop=FALSE])
  dropped <- meta2[[id_col]][!cc]; if (length(dropped)) meta2 <- meta2[cc, , drop=FALSE]
  mat2 <- mat[, meta2[[id_col]], drop=FALSE]
  X <- model.matrix(as.formula(paste("~", paste(design_vars, collapse=" + "))), data=meta2)
  list(mat=mat2, meta=meta2, X=X, dropped=dropped)
}

# ---------- main pipeline ----------
run_pipeline <- function(
  root_dir = ".",
  metadata_csv = "data/subject_metadata.csv",
  taxa_table = "data/metaphlan_rel_ab.csv",
  taxa_table_is_samples_in_columns = TRUE,
  taxa_feature_id_col = "clade_name",
  taxa_is_metaphlan = TRUE,
  taxa_select_level = "species",
  functional_table = "data/humann_pathway_counts.csv",
  functional_feature_id_col = "Pathway",
  enterotypes_csv = NULL,
  enterotypes_filter_true_only = FALSE,
  id_col = "sample_name",
  group_col = "Constipation",
  group_case_values = c("Y","Yes","TRUE"),
  diagnosis_col = "Case_status",
  diagnosis_pd_label = "PD",
  covariates = c("Age_at_collection","Sex","BMI","Laxatives","collection_method"),
  fig3_top_n = 20,
  fig3_make_volcano = TRUE,
  prevalence_threshold = 1e-4,
  prevalence_min_prop = 0.05,
  seed = 42,
  out_dir = "results"
) {
  set.seed(seed)
  o <- file.path(root_dir, out_dir); figs <- file.path(o, "figures"); tbls <- file.path(o, "tables"); supp <- file.path(tbls, "supplement")
  safe_mkdir(o); safe_mkdir(figs); safe_mkdir(tbls); safe_mkdir(supp)
  log_msg("Starting pipeline...")

  # metadata (PD-only)
  meta_raw <- fread(locate(root_dir, metadata_csv)) %>% as_tibble()
  stopifnot(id_col %in% names(meta_raw))
  if (!is.null(diagnosis_col) && diagnosis_col %in% names(meta_raw)) {
    meta_raw <- meta_raw %>% dplyr::filter(.data[[diagnosis_col]] == diagnosis_pd_label)
  }
  stopifnot(group_col %in% names(meta_raw))
  meta_raw$Group <- standardize_group_values(meta_raw[[group_col]], group_case_values)
  meta_raw <- meta_raw %>% dplyr::filter(!is.na(Group))

  # taxa table (MetaPhlAn-style RA)
  tax <- fread(locate(root_dir, taxa_table)) %>% as_tibble()
  stopifnot(taxa_feature_id_col %in% names(tax))
  if (isTRUE(taxa_is_metaphlan) && !is.null(taxa_select_level)) {
    clade <- tax[[taxa_feature_id_col]]
    keep <- switch(taxa_select_level,
      species = str_detect(clade, fixed("s__")),
      genus   = str_detect(clade, fixed("g__")) & !str_detect(clade, fixed("s__")),
      family  = str_detect(clade, fixed("f__")) & !str_detect(clade, fixed("g__")),
      order   = str_detect(clade, fixed("o__")) & !str_detect(clade, fixed("f__")),
      class   = str_detect(clade, fixed("c__")) & !str_detect(clade, fixed("o__")),
      phylum  = str_detect(clade, fixed("p__")) & !str_detect(clade, fixed("c__"))
    )
    tax <- tax[keep, , drop=FALSE]
  }
  counts0 <- if (isTRUE(taxa_table_is_samples_in_columns)) {
    tax %>% tibble::column_to_rownames(taxa_feature_id_col) %>% as.matrix()
  } else {
    t(tax %>% tibble::column_to_rownames(taxa_feature_id_col) %>% as.matrix())
  }

  # align ids
  ali <- align_counts_meta(counts0, meta_raw, id_col)
  counts <- ali$counts; meta <- ali$meta
  writeLines(ali$diagnostics, file.path(tbls, "id_alignment_report.txt"))
  if (ncol(counts) == 0) stop("No overlapping sample IDs after normalization. See tables/id_alignment_report.txt")

  # phyloseq object (taxonomy matrix must be character)
  OTU <- otu_table(counts, taxa_are_rows=TRUE)
  SAM <- sample_data(meta %>% column_to_rownames(id_col))
  TAX <- tax_table(matrix(NA_character_, nrow=nrow(counts), ncol=1, dimnames=list(rownames(counts), "feature")))
  ps <- phyloseq(OTU, SAM, TAX)
  ps_rel <- normalize_rel(ps)

  # --------- Table 1 ----------
  tab1 <- meta %>% count(Group) %>% bind_rows(tibble(Group="Total", n=nrow(meta))) %>% 
          select(Group, N=n) %>% arrange(match(Group, c("Total","control","case")))
  save_table_with_html(tab1, file.path(tbls, "Table1_sample_counts.csv"), file.path(tbls, "Table1_sample_counts.html"),
                       title="Table 1. PD-only cohort counts used in analyses",
                       footnotes=c("Definitions: control = non-constipated PD; case = constipated PD.",
                                   "Counts are after sample-ID alignment across metadata and abundance tables."))

  # --------- Fig 1: Alpha diversity ----------
  log_msg("Alpha diversity...")
  m_rel <- as(otu_table(ps_rel), "matrix"); m_rel[!is.finite(m_rel)] <- 0
  calc_shannon <- function(p) { p <- p[p>0]; -sum(p*log(p)) }
  calc_simpson <- function(p) 1 - sum(p^2)
  alpha_tbl <- tibble(!!id_col := colnames(m_rel),
                      Shannon = apply(m_rel, 2, calc_shannon),
                      Simpson = apply(m_rel, 2, calc_simpson),
                      Observed = colSums(m_rel > 0)) %>% left_join(meta, by=id_col)
  do_wilcox <- function(df, var) {
    x <- df %>% filter(Group=="case") %>% pull(.data[[var]])
    y <- df %>% filter(Group=="control") %>% pull(.data[[var]])
    wt <- wilcox.test(x, y, exact=FALSE)
    tibble(metric=var, p=wt$p.value, p_fmt=fmt_p(wt$p.value),
           med_case=median(x,na.rm=TRUE), med_ctrl=median(y,na.rm=TRUE))
  }
  alpha_stats <- bind_rows(do_wilcox(alpha_tbl,"Shannon"),
                           do_wilcox(alpha_tbl,"Simpson"),
                           do_wilcox(alpha_tbl,"Observed"))
  fwrite(alpha_stats, file.path(tbls, "alpha_diversity_stats.csv"))

  p_alpha_shannon  <- ggplot(alpha_tbl, aes(x=Group, y=Shannon, fill=Group)) +
    geom_boxplot(outlier.shape=NA, width=.6, color="black", linewidth=.3) +
    geom_jitter(color="black", width=.12, alpha=.6, size=1) +
    labs(x=NULL, y="Shannon diversity") + theme_base() + scale_fill_group()
  p_alpha_simpson  <- ggplot(alpha_tbl, aes(x=Group, y=Simpson, fill=Group)) +
    geom_boxplot(outlier.shape=NA, width=.6, color="black", linewidth=.3) +
    geom_jitter(color="black", width=.12, alpha=.6, size=1) +
    labs(x=NULL, y="Gini-Simpson (1 - D)") + theme_base() + scale_fill_group()
  p_alpha_observed <- ggplot(alpha_tbl, aes(x=Group, y=Observed, fill=Group)) +
    geom_boxplot(outlier.shape=NA, width=.6, color="black", linewidth=.3) +
    geom_jitter(color="black", width=.12, alpha=.6, size=1) +
    labs(x=NULL, y="Observed richness") + theme_base() + scale_fill_group()

  save_figure(p_alpha_shannon,  file.path(figs, "alpha_shannon"), width=3.2, height=4.0)
  save_figure(p_alpha_simpson,  file.path(figs, "alpha_simpson"), width=3.2, height=4.0)
  save_figure(p_alpha_observed, file.path(figs, "alpha_observed"), width=3.2, height=4.0)
  save_figure((p_alpha_shannon | p_alpha_simpson | p_alpha_observed) + plot_annotation(tag_levels='a'),
              file.path(figs, "Fig1_alpha"), width=9.6, height=4.0)

  # --------- Beta-diversity + PERMANOVA (Fig 2) ----------
  log_msg("Beta diversity + PERMANOVA + dispersion...")
  bray <- phyloseq::distance(ps_rel, method="bray")
  otu_raw <- as(otu_table(ps), "matrix"); clr_all <- t(clr_transform(otu_raw)); aitch <- dist(clr_all, method="euclidean")
  present_covs <- covariates[covariates %in% names(meta)]
  ad_bray  <- adonis2(bray  ~ ., data=meta %>% select(Group, all_of(present_covs)), permutations=999, by="terms")
  ad_aitch <- adonis2(aitch ~ ., data=meta %>% select(Group, all_of(present_covs)), permutations=999, by="terms")

  get_group_row <- function(ad) {
    df <- as.data.frame(ad); rn <- rownames(df)
    hit <- which(grepl("^Group$", rn)); if (!length(hit)) hit <- which(grepl("^Group", rn))
    if (!length(hit)) return(tibble(term="Group", R2=NA_real_, p=NA_real_))
    tibble(term="Group", R2=as.numeric(df[hit,"R2"]), p=as.numeric(df[hit,"Pr(>F)"]))
  }
  perma <- bind_rows(get_group_row(ad_bray)  %>% mutate(distance="Bray-Curtis"),
                     get_group_row(ad_aitch) %>% mutate(distance="Aitchison"))
  perma$q <- p.adjust(perma$p, method="BH"); perma$p_fmt <- fmt_p(perma$p); perma$q_fmt <- fmt_p(perma$q)
  fwrite(perma %>% select(distance, term, R2, p, p_fmt, q, q_fmt), file.path(tbls, "Table2_permanova_summary.csv"))

  bd_bray  <- betadisper(bray, meta$Group);  disp_bray  <- permutest(bd_bray, permutations=999)
  bd_aitch <- betadisper(aitch, meta$Group); disp_aitch <- permutest(bd_aitch, permutations=999)
  capture.output(list(bray=disp_bray, aitchison=disp_aitch), file=file.path(tbls, "permdisp_tests.txt"))

  ord_plot <- function(d, md) {
    ord <- cmdscale(d, k=2, eig=TRUE); eig <- ord$eig / sum(ord$eig); eig[!is.finite(eig)] <- 0
    df <- as.data.frame(ord$points); names(df) <- c("Axis1","Axis2"); df$sample_id <- rownames(df)
    df <- df %>% left_join(md %>% rename(sample_id=!!id_col), by="sample_id")
    ggplot(df, aes(Axis1, Axis2, colour=Group)) +
      geom_point(alpha=.85, size=1.5) +
      labs(x=paste0("Axis 1 (", percent(pmax(eig[1],0)), ")"),
           y=paste0("Axis 2 (", percent(pmax(eig[2],0)), ")")) +
      theme_base() + scale_colour_group()
  }
  p_pcoa_bray  <- ord_plot(bray,  meta)
  p_pcoa_aitch <- ord_plot(aitch, meta)
  cap_txt <- sprintf("PERMANOVA R^2 - Bray: %s (q %s) | Aitchison: %s (q %s)",
                     {x<-perma %>% filter(distance=="Bray-Curtis") %>% pull(R2); ifelse(length(x), sprintf("%.3f",x), "NA")},
                     {x<-perma %>% filter(distance=="Bray-Curtis") %>% pull(q_fmt); ifelse(length(x), x, "NA")},
                     {x<-perma %>% filter(distance=="Aitchison") %>% pull(R2); ifelse(length(x), sprintf("%.3f",x), "NA")},
                     {x<-perma %>% filter(distance=="Aitchison") %>% pull(q_fmt); ifelse(length(x), x, "NA")})
  save_figure(p_pcoa_bray,  file.path(figs, "pcoa_bray"),      width=3.6, height=3.6)
  save_figure(p_pcoa_aitch, file.path(figs, "pcoa_aitchison"), width=3.6, height=3.6)
  save_figure((p_pcoa_bray | p_pcoa_aitch) + plot_annotation(tag_levels='a', caption=cap_txt),
              file.path(figs, "Fig2_beta"), width=7.2, height=3.5)

  # --------- Fig 3: Differential taxa (CLR + limma) ----------
  log_msg("Differential taxa (CLR + limma)...")
  taxa_mat <- as(otu_table(ps), "matrix")
  rel_m0 <- sweep(taxa_mat, 2, colSums(taxa_mat), "/"); rel_m0[!is.finite(rel_m0)] <- 0
  keep_prev <- rowMeans(rel_m0 > prevalence_threshold, na.rm=TRUE) >= prevalence_min_prop
  if (any(keep_prev)) {
    taxa_mat_f <- taxa_mat[keep_prev, , drop=FALSE]
    clr_m <- clr_transform(taxa_mat_f)
    keep_var <- apply(clr_m, 1, function(r) sd(r, na.rm=TRUE) > 0)
    clr_m <- clr_m[keep_var, , drop=FALSE]
    prep <- prepare_design(meta, group_col="Group", covars=covariates); meta_d <- prep$meta; covs <- prep$covars
    aa <- align_for_model(clr_m, meta_d, id_col, covs)
    if (length(aa$dropped)) writeLines(c("Dropped for NA in design vars:", paste0("  - ", aa$dropped)), file.path(tbls, "design_drop_report.txt"))
    fit <- suppressWarnings(lmFit(aa$mat, aa$X))  # suppress harmless "Partial NA coefficients"
    ok_rows <- rowSums(is.na(fit$coefficients)) == 0
    if (!all(ok_rows)) { fit <- fit[ok_rows, ]; aa$mat <- aa$mat[ok_rows, , drop=FALSE] }
    fit <- eBayes(fit)
    grp_col <- grep("^Group", colnames(aa$X)); stopifnot(length(grp_col) == 1)
    tt <- topTable(fit, coef=grp_col, number=Inf, sort.by="P") %>% rownames_to_column("feature") %>% as_tibble() %>%
      mutate(adjP_fmt=fmt_p(adj.P.Val), P_fmt=fmt_p(P.Value)) %>% filter(is.finite(logFC))
    fwrite(tt, file.path(tbls, "taxa_differential_limma_clr.csv"))
    ntax_tested <- nrow(tt)

    tt <- tt %>% mutate(feature_short=short_tax_label(feature),
                        dir=factor(ifelse(logFC>=0,"case","control"), levels=c("control","case")))
    top20 <- tt %>% arrange(adj.P.Val, desc(abs(logFC))) %>% head(fig3_top_n)

    p_vol <- ggplot(tt %>% filter(!is.na(adj.P.Val)), aes(logFC, -log10(adj.P.Val), colour=dir)) +
      geom_point(alpha=.7, size=1.2) +
      geom_vline(xintercept=c(-0.1,0.1), linetype=2, linewidth=.3) +
      geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=.3) +
      labs(x="logFC (case/control)", y="-log10(q)") + theme_base() + scale_colour_group() + guides(colour="none")

    p_top <- top20 %>%
      mutate(feature_short=fct_reorder(feature_short, logFC)) %>%
      mutate(feature_short=fct_relabel(feature_short, ~ifelse(.x=="Flavonifractor_plauti","Flavonifractor_plautii",.x))) %>%
      ggplot(aes(x=feature_short, y=logFC, fill=dir)) +
      geom_col(width=0.7, color="black", linewidth=0.15) +
      coord_flip() +
      labs(x=NULL, y="Effect size (logFC, CLR)") +
      theme_base() + scale_fill_group() + guides(fill="none") +
      theme(axis.text.y=element_text(size=base_font_size-1.5), plot.margin=margin(10,18,10,10))

    save_figure((p_vol | p_top) + plot_annotation(tag_levels='a'), file.path(figs, "Fig3_taxa"), width=7.6, height=4.2)

    prev <- rowMeans(rel_m0[rownames(aa$mat), , drop=FALSE] > prevalence_threshold)
    mean_control <- rowMeans(rel_m0[rownames(aa$mat), aa$meta$Group=="control", drop=FALSE], na.rm=TRUE)
    mean_case    <- rowMeans(rel_m0[rownames(aa$mat), aa$meta$Group=="case",    drop=FALSE], na.rm=TRUE)
    tt2 <- tt %>% mutate(Prevalence=signif(prev[feature],3),
                         Mean_RA_control=signif(mean_control[feature],3),
                         Mean_RA_case=signif(mean_case[feature],3))
    top20_tbl <- tt2 %>% arrange(adj.P.Val, desc(abs(logFC))) %>% transmute(
      Rank=row_number(), Taxon=feature_short, logFC=signif(logFC,4), q=signif(adj.P.Val,3),
      Direction=ifelse(logFC>=0, "up case", "up control"),
      Prevalence=signif(Prevalence,3), Mean_RA_control=signif(Mean_RA_control,3), Mean_RA_case=signif(Mean_RA_case,3)
    ) %>% head(fig3_top_n)
    t2_foot <- c(sprintf("Benjamini-Hochberg FDR; features tested: %d (after prevalence filter RA>%g in >=%g%% of samples).",
                         ntax_tested, prevalence_threshold, 100*prevalence_min_prop),
                 sprintf("Covariates: %s.", paste(covs, collapse=", ")),
                 "Sign: positive logFC = higher in case (constipated PD).")
    save_table_with_html(top20_tbl, file.path(tbls, "Table2_top20_taxa.csv"), file.path(tbls, "Table2_top20_taxa.html"),
                         title="Table 2. Top 20 differentially abundant taxa (case vs control, PD-only)", footnotes=t2_foot)
  } else {
    writeLines("No taxa pass prevalence filter.", file.path(tbls, "taxa_differential_error.txt"))
  }

  # --------- Fig 4: Pathways (voom-limma) ----------
  if (!is.null(functional_table)) {
    if (!file.exists(locate(root_dir, functional_table))) {
      log_msg("Functional table not found; skipping pathways.")
    } else {
      log_msg("Differential pathways (voom-limma)...")
      func <- fread(locate(root_dir, functional_table)) %>% as_tibble()
      stopifnot(functional_feature_id_col %in% names(func))
      mat0 <- func %>% tibble::column_to_rownames(functional_feature_id_col) %>% as.matrix()
      keep_samp <- intersect(colnames(mat0), meta[[id_col]])
      mat <- mat0[, keep_samp, drop=FALSE]
      meta_f <- meta %>% filter(.data[[id_col]] %in% keep_samp) %>% arrange(match(.data[[id_col]], keep_samp))
      prep_f <- prepare_design(meta_f, group_col="Group", covars=covariates); meta_f <- prep_f$meta; covs_f <- prep_f$covars
      aa <- align_for_model(mat, meta_f, id_col, covs_f)
      if (length(aa$dropped)) writeLines(c("Pathways: dropped for NA in design vars:", paste0("  - ", aa$dropped)), file.path(tbls, "pathways_design_drop_report.txt"))

      y <- DGEList(counts=round(aa$mat))
      v <- voom(y, design=aa$X, plot=FALSE)
      fit <- suppressWarnings(lmFit(v, aa$X))
      ok_rows <- rowSums(is.na(fit$coefficients)) == 0
      if (!all(ok_rows)) { fit <- fit[ok_rows, ]; v$E <- v$E[ok_rows, , drop=FALSE] }
      fit <- eBayes(fit)
      grp_col <- grep("^Group", colnames(aa$X)); stopifnot(length(grp_col) == 1)
      tt <- topTable(fit, coef=grp_col, number=Inf, sort.by="P") %>% rownames_to_column("feature") %>% as_tibble() %>%
        mutate(adjP_fmt=fmt_p(adj.P.Val), P_fmt=fmt_p(P.Value))
      fwrite(tt, file.path(tbls, "pathways_differential_limma.csv"))
      n_pwy_tested <- nrow(tt)

      tt <- tt %>% mutate(feature_short=short_path_label(feature),
                          dir=factor(ifelse(logFC>=0,"case","control"), levels=c("control","case")))
      pwy_top <- tt %>% arrange(adj.P.Val, desc(abs(logFC))) %>% head(20) %>% mutate(feature_short=fct_reorder(feature_short, logFC))
      xmax <- max(abs(pwy_top$logFC), na.rm=TRUE); lims <- c(-xmax, xmax) * 1.1
      p_pwy <- ggplot(pwy_top, aes(x=feature_short, y=logFC, fill=dir)) +
        geom_col(width=0.7, color="black", linewidth=0.15) +
        coord_flip() +
        scale_y_continuous(limits=lims, breaks=breaks_extended(5), expand=expansion(mult=c(0.02,0.05))) +
        labs(x=NULL, y="Effect size (log2FC)") +
        theme_base() + scale_fill_group() + guides(fill="none") +
        theme(axis.text.y=element_text(size=base_font_size-1.5), plot.margin=margin(10,18,10,10))
      save_figure(p_pwy + plot_annotation(tag_levels='a'), file.path(figs, "Fig4_pathways"), width=8.8, height=5.0)

      nf <- calcNormFactors(y); cpm_mat <- cpm(y, normalized.lib.sizes=TRUE, log=FALSE)
      mean_cpm_control <- rowMeans(cpm_mat[, aa$meta$Group=="control", drop=FALSE], na.rm=TRUE)
      mean_cpm_case    <- rowMeans(cpm_mat[, aa$meta$Group=="case",    drop=FALSE], na.rm=TRUE)
      tbl3 <- tt %>% arrange(adj.P.Val, desc(abs(logFC))) %>% head(20) %>% transmute(
        Rank=row_number(), Pathway_ID=feature, Description=feature_short, logFC=signif(logFC,4), q=signif(adj.P.Val,3),
        Direction=ifelse(logFC>=0,"up case","up control"),
        Mean_CPM_control=signif(mean_cpm_control[feature],4), Mean_CPM_case=signif(mean_cpm_case[feature],4)
      )
      t3_foot <- c(sprintf("Benjamini-Hochberg FDR; pathways tested: %d.", n_pwy_tested),
                   sprintf("Covariates: %s.", paste(covs_f, collapse=", ")),
                   "Sign: positive log2FC = higher in case (constipated PD).",
                   "Pathway table derived from published HUMAnN outputs (not re-run).")
      save_table_with_html(tbl3, file.path(tbls, "Table3_top20_pathways.csv"), file.path(tbls, "Table3_top20_pathways.html"),
                           title="Table 3. Top 20 differentially abundant pathways (case vs control, PD-only)", footnotes=t3_foot)
    }
  }

  # --------- Fig 5: Enterotypes ----------
  log_msg("Enterotypes (CLR -> PCA -> clustering)...")
  otu_raw <- as(otu_table(ps), "matrix"); clr <- t(clr_transform(otu_raw)); clr[!is.finite(as.matrix(clr))] <- 0
  pca <- prcomp(clr, center=TRUE, scale.=FALSE); var_expl <- cumsum(pca$sdev^2)/sum(pca$sdev^2)
  n_pc <- min(which(var_expl >= 0.8)[1], 10); if (is.na(n_pc) || n_pc < 2) n_pc <- min(10, ncol(pca$x))
  pcs <- pca$x[, 1:n_pc, drop=FALSE]
  clust_ok <- TRUE; ent_df <- NULL; sel_tbl <- NULL; clust_msg <- NULL
  tryCatch({
    mc <- Mclust(pcs, G=1:6, verbose=FALSE)
    z <- mc$classification
    # Robust BIC tidy (handle matrix or vector; never error)
    bic <- mc$BIC
    sel_tbl <<- NULL
    if (!is.null(bic)) {
      if (is.matrix(bic)) {
        df <- as.data.frame(bic); df$k <- suppressWarnings(as.integer(rownames(df)))
        sel_tbl <<- tidyr::pivot_longer(df, cols = -k, names_to="model", values_to="BIC")
      } else {
        kv <- suppressWarnings(as.integer(names(bic))); if (any(is.na(kv))) kv <- seq_along(bic)
        sel_tbl <<- tibble(k = kv, model = NA_character_, BIC = as.numeric(bic))
      }
    }
    ent_df <<- tibble(!!id_col := rownames(pcs), enterotype=factor(z, levels=sort(unique(z)))) %>% left_join(meta, by=id_col)
  }, error=function(e){ clust_ok<<-FALSE; clust_msg<<-conditionMessage(e) })
  if (!clust_ok || is.null(ent_df)) {
    log_msg(paste("mclust failed; fallback to kmeans:", clust_msg))
    ks <- 2:6; sil_vals <- rep(NA_real_, length(ks)); set.seed(seed)
    for (i in seq_along(ks)) {
      km <- try(kmeans(pcs, centers=ks[i], nstart=25), silent=TRUE); if (inherits(km,"try-error")) next
      sil <- silhouette(km$cluster, dist(pcs)); sil_vals[i] <- mean(sil[,3])
    }
    if (!all(is.na(sil_vals))) {
      k_best <- ks[which.max(sil_vals)]; km <- kmeans(pcs, centers=k_best, nstart=50)
      sel_tbl <- tibble(k=ks, metric="silhouette", value=sil_vals)
      ent_df  <- tibble(!!id_col := rownames(pcs), enterotype=factor(km$cluster, levels=sort(unique(km$cluster)))) %>% left_join(meta, by=id_col)
    } else {
      writeLines("Clustering failed (mclust and kmeans).", file.path(tbls, "enterotype_error.txt"))
    }
  }
  if (!is.null(ent_df)) {
    tab <- ent_df %>% count(Group, enterotype) %>% tidyr::complete(Group, enterotype, fill=list(n=0))
    p_fisher <- tryCatch(fisher.test(table(ent_df$Group, ent_df$enterotype))$p.value, error=function(e) NA_real_)
    fwrite(tab %>% pivot_wider(names_from=Group, values_from=n), file.path(tbls, "Table5_enterotype_counts.csv"))
    writeLines(sprintf("Fisher exact test P = %s", fmt_p(p_fisher)), file.path(tbls, "enterotype_fisher.txt"))
    p_ent <- ggplot(tab, aes(x=enterotype, y=n, fill=Group)) + geom_col(position="dodge") +
      labs(x=NULL, y="Sample count") + theme_base() + scale_fill_group()
    save_figure(p_ent + plot_annotation(tag_levels='a'), file.path(figs, "Fig5_enterotypes"), width=7.2, height=3.8)
    if (!is.null(sel_tbl)) fwrite(sel_tbl, file.path(tbls, "enterotype_model_selection.csv"))
    fwrite(ent_df, file.path(tbls, "enterotype_assignments.csv"))
  }

  # --------- Supplement + session info ----------
  fwrite(alpha_stats, file.path(supp, "S1_alpha_diversity_stats.csv"))
  if (file.exists(file.path(tbls, "taxa_differential_limma_clr.csv"))) fwrite(fread(file.path(tbls, "taxa_differential_limma_clr.csv")), file.path(supp, "S2_taxa_full.csv"))
  if (file.exists(file.path(tbls, "pathways_differential_limma.csv"))) fwrite(fread(file.path(tbls, "pathways_differential_limma.csv")), file.path(supp, "S3_pathways_full.csv"))
  fwrite(perma, file.path(supp, "S4_permanova_summary.csv"))
  if (exists("ent_df") && !is.null(ent_df)) fwrite(ent_df, file.path(supp, "S5_enterotype_assignments.csv"))
  fwrite(tibble(id_alignment_report = ali$diagnostics), file.path(supp, "S6_id_alignment_report.csv"))
  writeLines(capture.output(sessionInfo()), file.path(tbls, "session_info.txt"))
  writeLines(c(
    "Main Figures: Fig1_alpha, Fig2_beta, Fig3_taxa, Fig4_pathways, Fig5_enterotypes (PNG + PDF).",
    "Main Tables: Table1_sample_counts (CSV+HTML), Table2_top20_taxa (CSV+HTML), Table3_top20_pathways (CSV+HTML).",
    "Alpha metrics computed on relative-abundance; Observed = count of taxa with RA>0.",
    "PERMANOVA by terms with covariates; dispersion via betadisper.",
    "Taxa: CLR + limma with prevalence filter; Pathways: voom-limma.",
    "Covariates: Age_at_collection, Sex, BMI, Laxatives, collection_method.",
    "Supplement in tables/supplement (S1-S6)."
  ), file.path(o, "README_outputs.txt"))

  log_msg("Pipeline complete. Outputs in:", normalizePath(o))
  invisible(TRUE)
}
