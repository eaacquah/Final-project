library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Dr.eg.db)

# Setting paths and parameters 
base <- "C:/Users/USER/Downloads"
setwd(base)
samples <- c("SRR24874370", "SRR24874372", "SRR24874373")
genotypes <- c("WT", "P23H", "P23H")
celltypes_to_test <- c("Rods", "Cones", "Microglia", "Muller_glia", "RGC",
                       "Progenitors", "Amacrine", "Bipolar", "Horizontal", "RPE")

seurat_list <- lapply(seq_along(samples), function(i) {
  
  sample <- samples[i]
  
  counts <- ReadMtx(
    mtx = file.path(base, sample, "matrix.mtx"),
    features = file.path(base, sample, "features.tsv"),
    cells = file.path(base, sample, "barcodes.tsv"),
    feature.column = 1
  )
  
  obj <- CreateSeuratObject(
    counts,
    project = sample,
    min.cells = 3,
    min.features = 200
  )
  
  obj$genotype <- genotypes[i]
  obj$sample <- sample
  
  obj
})

names(seurat_list) <- samples

combined <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = samples,
  project = "P23H_Retina"
)

saveRDS(combined, "P23H_Seurat_raw.rds")

# Inspecting the combined Seurat object
combined
dim(combined)

# Checking sample identities
table(combined$sample, combined$genotype)

# Calculating mitochondrial gene percentage
# a. Searching for mitochondrial genes
mt_genes <- grep("^mt[-.]", rownames(combined), value = TRUE, ignore.case = TRUE)
if (!length(mt_genes)) stop("No mitochondrial gene symbols found. Check features.tsv.")

# b. Filtering mitochondrial genes
combined[["percent.mt"]] <- PercentageFeatureSet(combined, features = mt_genes)
summary(combined$percent.mt)

# Visualizing QC distributions
# a. Violin plots by sample
p_plot <- p_plot <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            group.by = "sample", ncol = 3, pt.size = 0)
print(p_plot)
ggsave("QC_Violin_After_Filtering.pdf", p_plot, width = 14, height = 5)
ggsave("QC_Violin_After_Filtering.tiff", p_plot, width = 14, height = 5, dpi = 300, compression = "lzw")
print(p_plot)
ggsave("QC_Violin_Before_Filtering.pdf", p_plot, width = 14, height = 5)
ggsave("QC_Violin_Before_Filtering.tiff", p_plot, width = 14, height = 5, dpi = 300, compression = "lzw")

# b. Scatter plots
p_plot <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p_plot)
ggsave("QC_Scatter_Count_vs_Features.pdf", p_plot, width = 7, height = 5)
ggsave("QC_Scatter_Count_vs_Features.tiff", p_plot, width = 7, height = 5, dpi = 300, compression = "lzw")
p_plot <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(p_plot)
ggsave("QC_Scatter_Count_vs_Mito.pdf", p_plot, width = 7, height = 5)
ggsave("QC_Scatter_Count_vs_Mito.tiff", p_plot, width = 7, height = 5, dpi = 300, compression = "lzw")

# Summary statistics
qc_summary <- combined@meta.data %>%
  group_by(sample, genotype) %>%
  summarise(n_cells = n(), median_genes = median(nFeature_RNA),
            median_UMIs = median(nCount_RNA), median_mito = median(percent.mt),
            .groups = "drop")
print(qc_summary)
write.csv(qc_summary, "QC_summary_before_filtering.csv", row.names = FALSE)

# Filtering low quality cells
# a. Recording cells before filtering
cells_before <- table(combined$sample)

# b. QC filtering
combined <- subset(combined, subset = nFeature_RNA > 300 &
                     nFeature_RNA < 6000 & percent.mt < 10)

# Checking retained cells
print(cells_before)
print(table(combined$sample))

# Visualizing after filtering
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", ncol = 3, pt.size = 0)

# Saving filtered object
saveRDS(combined, "combined_postQC.rds")

# SCTransform normalization
combined <- SCTransform(combined, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(combined, "combined_SCT.rds")

# PCA, UMAP and clustering
# PCA
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)

# a. Determining suitable PC range
p_plot <- ElbowPlot(combined, ndims = 50)
print(p_plot)
ggsave("PCA_ElbowPlot.pdf", p_plot, width = 7, height = 5)
ggsave("PCA_ElbowPlot.tiff", p_plot, width = 7, height = 5, dpi = 300, compression = "lzw")

# b. Initial choice of 25 PCs
combined <- RunUMAP(combined, dims = 1:25)
combined <- FindNeighbors(combined, dims = 1:25)
combined <- FindClusters(combined, resolution = 1.0)

# Inspecting cluster distribution
table(combined$seurat_clusters)

# Save
saveRDS(combined, "combined_after_DR.rds")

# Visualizing UMAPs
# a. Per clusters
p_plot <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
print(p_plot)
ggsave("UMAP_Clusters.pdf", p_plot, width = 9, height = 7)
ggsave("UMAP_Clusters.tiff", p_plot, width = 9, height = 7, dpi = 300, compression = "lzw")

# b. Per samples
p_plot <- DimPlot(combined, reduction = "umap", group.by = "sample")
print(p_plot)
ggsave("UMAP_Samples.pdf", p_plot, width = 9, height = 7)
ggsave("UMAP_Samples.tiff", p_plot, width = 9, height = 7, dpi = 300, compression = "lzw")

# c. Per genotype overlap
p_plot <- DimPlot(combined, reduction = "umap", group.by = "genotype", shuffle = TRUE)
print(p_plot)
ggsave("UMAP_Genotype.pdf", p_plot, width = 9, height = 7)
ggsave("UMAP_Genotype.tiff", p_plot, width = 9, height = 7, dpi = 300, compression = "lzw")

# Genotypes in separate panels
p_plot <- DimPlot(combined, reduction = "umap", split.by = "genotype")
print(p_plot)
ggsave("UMAP_Genotype_Split.pdf", p_plot, width = 12, height = 6)
ggsave("UMAP_Genotype_Split.tiff", p_plot, width = 12, height = 6, dpi = 300, compression = "lzw")

# Identifying cluster marker genes
combined <- JoinLayers(combined, assay = "RNA")
combined <- NormalizeData(combined, assay = "RNA", verbose = FALSE)
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "seurat_clusters"

# Identifying positive markers
cluster_markers <- FindAllMarkers(combined, assay = "RNA", only.pos = TRUE,
                                  min.pct = 0.1, logfc.threshold = 0.25,
                                  test.use = "wilcox")

# Save results
write.csv(cluster_markers, "cluster_marker_genes.csv", row.names = FALSE)

# Top five markers per cluster
top5 <- cluster_markers %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>% ungroup()
print(top5)
write.csv(top5, "top5_cluster_markers.csv", row.names = FALSE)

# Annotating retinal cell types
retina_markers <- list(
  Amacrine = c("cabp1a", "calb2a", "tfap2a", "tfap2b", "syt2a", "chata", "gad1b", "gad2", "slc6a9a", "ptf1a", "foxn4", "neurod6b"),
  Bipolar = c("cabp2a", "prkca", "vsx1", "vsx2", "grm6a", "grm6b", "rs1a", "bhlhe23", "pcp2a", "otx2b", "isl1"),
  Rods = c("gnat1", "rho", "pde6a", "pde6b", "rom1a", "rom1b", "sagaa", "sagab", "nr2e3"),
  Cones = c("gnat2", "pde6c", "arr3a", "arr3b", "opn1lw1", "opn1mw1", "opn1sw1", "opn1sw2", "clul1", "crx", "thrb"),
  Horizontal = c("prox1a", "lim1a", "calb1a", "cx52.6", "cx55.5", "cx52.9", "onecut1", "onecut2", "lhx1a"),
  Progenitors = c("her4.3", "insm1a", "her15.1", "ascl1a", "pcna", "mki67", "cdkn1ca", "pax6a", "ccnd1"),
  RPE = c("rpe65a", "dct", "pmela", "pmelb", "tyr", "tyrp1b", "mitfa", "best1", "rlbp1b", "slc16a8a"),
  RGC = c("rbpms2a", "rbpms2b", "pou4f1", "isl2b", "pou4f2", "pou4f3", "atoh7", "sncg", "robo2", "gap43"),
  Oligodendrocyte = c("plp1b", "cd59", "cd9b", "mbpa", "mpz", "sox10", "olig2"),
  Muller_glia = c("gfap", "cahz", "ptgdsb.1", "glula", "slc1a3a", "rlbp1a", "vim", "fabp7a", "slc1a2b", "sox2", "aldh1a3"),
  Microglia = c("cd74a", "cd74b", "apoc1", "mpeg1.1", "c1qa", "c1qb", "csf1ra", "p2ry12", "hexb", "ptprc", "slc7a7"),
  Erythroid = c("hbaa2", "hbba2", "slc4a1a", "hemgn")
)

# Checking marker availability
# a. Keeping markers present in the expression matrix
retina_markers_clean <- lapply(retina_markers, intersect, y = rownames(combined))
# b. Removing empty groups
retina_markers_clean <- Filter(length, retina_markers_clean)
if (!length(retina_markers_clean)) stop("No retinal marker symbols matched the matrix.")

# Dot plot
p_plot <- DotPlot(combined, features = retina_markers_clean, assay = "RNA") + RotatedAxis()
print(p_plot)
ggsave("Retinal_Marker_DotPlot.pdf", p_plot, width = 16, height = 8)
ggsave("Retinal_Marker_DotPlot.tiff", p_plot, width = 16, height = 8, dpi = 300, compression = "lzw")

# Assigning annotations
annotation_file <- "cluster_annotations.csv"
if (!file.exists(annotation_file)) {
  cluster_annotations <- data.frame(cluster = levels(combined$seurat_clusters),
                                    celltype = "Unassigned")
  write.csv(cluster_annotations, annotation_file, row.names = FALSE)
  stop("Fill in celltype values in cluster_annotations.csv, then rerun from 'Reading manually assigned annotations'.")
}

# Reading manually assigned annotations
cluster_annotations <- read.csv(annotation_file, stringsAsFactors = FALSE)
annotation_map <- setNames(cluster_annotations$celltype, as.character(cluster_annotations$cluster))
combined$celltype <- unname(annotation_map[as.character(combined$seurat_clusters)])

# Plotiing annotated UMAP
p_plot <- DimPlot(combined, group.by = "celltype", label = TRUE, repel = TRUE)
print(p_plot)
ggsave("UMAP_Annotated_Celltypes.pdf", p_plot, width = 10, height = 8)
ggsave("UMAP_Annotated_Celltypes.tiff", p_plot, width = 10, height = 8, dpi = 300, compression = "lzw")
saveRDS(combined, "annotated_retina_seurat.rds")

# Comparing cell-type proportions
cell_props <- combined@meta.data %>%
  count(sample, genotype, celltype, name = "count") %>%
  group_by(sample) %>% mutate(percent = 100 * count / sum(count)) %>% ungroup()

# Plotting cell-type proportions
p_plot <- ggplot(cell_props, aes(x = sample, y = percent, fill = celltype)) +
  geom_col() + theme_classic() +
  labs(x = "Sample", y = "Percentage of retained cells", fill = "Cell type")
print(p_plot)
ggsave("Celltype_Proportions.pdf", p_plot, width = 9, height = 6)
ggsave("Celltype_Proportions.tiff", p_plot, width = 9, height = 6, dpi = 300, compression = "lzw")
write.csv(cell_props, "celltype_proportions.csv", row.names = FALSE)


# 1. OVERALL DIFFERENTIAL GENE EXPRESSION: P23H vs WT 

Idents(combined) <- "genotype"
overall_DEGs <- FindMarkers(combined, ident.1 = "P23H", ident.2 = "WT",
                            assay = "RNA", test.use = "wilcox",
                            min.pct = 0.1, logfc.threshold = 0.25)
write.csv(overall_DEGs, "overall_P23H_vs_WT_exploratory_DEGs.csv")

# Reusable DEG processing and plots
prepare_de <- function(de) {
  de$gene <- rownames(de)
  de %>% filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
    mutate(regulation = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE ~ "Not significant"),
      neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)))
}
select_top <- function(de) {
  bind_rows(
    de %>% filter(regulation == "Upregulated") %>%
      slice_max(avg_log2FC, n = 10, with_ties = FALSE),
    de %>% filter(regulation == "Downregulated") %>%
      slice_min(avg_log2FC, n = 10, with_ties = FALSE))
}
make_volcano <- function(de, top, title, subtitle) {
  ggplot(de, aes(avg_log2FC, neg_log10_padj, color = regulation)) +
    geom_point(alpha = 0.65, size = 1.5) +
    geom_text_repel(data = top, aes(label = gene), size = 3.5,
                    color = "black", max.overlaps = Inf, seed = 123,
                    show.legend = FALSE) +
    scale_color_manual(values = c(Upregulated = "#D62728",
                                  Downregulated = "#1F77B4",
                                  `Not significant` = "gray70")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    theme_classic(base_size = 13) +
    labs(title = title, subtitle = subtitle,
         x = "Log2 Fold Change (P23H / WT)",
         y = "-Log10 Adjusted p-value", color = "Expression") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
}
make_deg_dot <- function(obj, genes, title) {
  DotPlot(obj, features = genes, group.by = "genotype", assay = "RNA", scale = TRUE) +
    RotatedAxis() + theme_classic(base_size = 12) +
    labs(title = title, x = "Genes", y = "Genotype") +
    theme(plot.title = element_text(hjust = 0.5))
}

# Plotting top 10 overall differentially expressed genes
# a. Adding gene names
volcano_df <- prepare_de(overall_DEGs)

# b. Filtering DEGs: Top 10 upregulated genes in P23H and top 10 downregulated genes in P23H
# Combine gene lists
top20 <- select_top(volcano_df)

# c. Displaying selected genes
print(top20)

# d. Generating and saving dot plot
if (nrow(top20)) {
  p_dot <- make_deg_dot(combined, top20$gene, "Top 10 Upregulated and Downregulated Genes")
  print(p_dot)
  
  ggsave("Top20_Directional_DEGs_DotPlot.pdf", p_dot, width = 15, height = 5)
  ggsave("Top20_Directional_DEGs_DotPlot.tiff", p_dot, width = 15, height = 5, dpi = 300, compression = "lzw")
}

# Volcano plot of overall differential gene expression: P23H vs WT zebrafish retina
# Labeling top 10 upregulated and top 10 downregulated genes
# a. Preparing volcano plot data: Counting genes in each category
print(table(volcano_df$regulation))

# b. Selecting top 10 upregulated genes and top 10 downregulated genes
# Combining genes for labeling and selected genes
print(top20 %>% dplyr::select(gene, avg_log2FC, p_val_adj, regulation))

# c. Creating labeled volcano plot
p_volcano <- make_volcano(volcano_df, top20,
                          "P23H vs WT Retinal Differential Expression",
                          "Exploratory tissue-level analysis")
# d. Displaying and saving volcano plot
print(p_volcano)

ggsave("Overall_P23H_vs_WT_Labeled_VolcanoPlot.pdf", p_volcano,
       width = 11, height = 8)
ggsave("Overall_P23H_vs_WT_Labeled_VolcanoPlot.tiff", p_volcano,
       width = 11, height = 8, dpi = 300, compression = "lzw")


# 2. CELL TYPE-SPECIFIC DIFFERENTIAL GENE EXPRESSION: P23H vs WT

# Creating output directories
for (folder in c("Celltype_DEG_Plots/DotPlots", "Celltype_DEG_Plots/VolcanoPlots",
                 "Celltype_DEG_Plots/Top20_Gene_Tables"))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

# Looping through each cell type
cell_degs <- list()
for (ct in celltypes_to_test) {
  if (!ct %in% combined$celltype) next
  
  # a. Creating subset seurat objects
  cells_ct <- subset(combined, subset = celltype == ct)
  Idents(cells_ct) <- "genotype"
  group_counts <- table(Idents(cells_ct))
  if (!all(c("WT", "P23H") %in% names(group_counts)) ||
      any(group_counts[c("WT", "P23H")] < 3)) {
    message("Skipping ", ct, ": insufficient cells or missing genotype")
    next
  }
  de_raw <- FindMarkers(cells_ct, ident.1 = "P23H", ident.2 = "WT",
                        assay = "RNA", test.use = "wilcox",
                        min.pct = 0.1, logfc.threshold = 0.25)
  cell_degs[[ct]] <- de_raw
  write.csv(de_raw, paste0(ct, "_P23H_vs_WT_exploratory_DEGs.csv"))
  de <- prepare_de(de_raw)
  if (!nrow(de)) next
  
  # b. Filtering DEGs: Top 10 upregulated in P23H and top 10 downregulated in P23H
  # Combining selected genes
  top20_table <- select_top(de)
  
  # c. Saving top 20 DEG table 
  write.csv(top20_table,
            file.path("Celltype_DEG_Plots/Top20_Gene_Tables", paste0(ct, "_Top20_DEGs.csv")),
            row.names = FALSE)
  # Ensuring selected genes exist
  top20_genes <- intersect(top20_table$gene, rownames(cells_ct))
  
  # d. Generating dot plot
  if (length(top20_genes)) {
    p_dot <- make_deg_dot(cells_ct, top20_genes,
                          paste0(ct, ": Top Upregulated and Downregulated Genes"))
    # Display in RStudio
    print(p_dot)
    # Save PDF
    ggsave(file.path("Celltype_DEG_Plots/DotPlots",
                     paste0(ct, "_Top20_DEGs_DotPlot.pdf")),
           p_dot, width = 15, height = 5)
    ggsave(file.path("Celltype_DEG_Plots/DotPlots",
                     paste0(ct, "_Top20_DEGs_DotPlot.tiff")),
           p_dot, width = 15, height = 5, dpi = 300, compression = "lzw")
  }
  
  # e. Preparing volcano plot data: Count upregulated and downregulated genes
  n_up <- sum(de$regulation == "Upregulated")
  n_down <- sum(de$regulation == "Downregulated")
  
  # f. Generating labeled volcano plot 
  # Selecting genes for volcano labeling; Labeling and extracting coordinates for the same top 20 genes used in the dot plot
  p_volcano <- make_volcano(de, top20_table, paste0(ct, ": P23H vs WT"),
                            paste0("Upregulated: ", n_up, " | Downregulated: ", n_down))
  
  # g. Displaying and saving volcano plots
  print(p_volcano)
  
  ggsave(file.path("Celltype_DEG_Plots/VolcanoPlots",
                   paste0(ct, "_Labeled_VolcanoPlot.pdf")),
         p_volcano, width = 11, height = 8)
  ggsave(file.path("Celltype_DEG_Plots/VolcanoPlots",
                   paste0(ct, "_Labeled_VolcanoPlot.tiff")),
         p_volcano, width = 11, height = 8, dpi = 300, compression = "lzw")
}

# Notice of completion
message("All available cell-type DEG plots completed.")

# Saving annotated retinal dataset
saveRDS(combined, "annotated_retina_final.rds")



# 3. ZEBRAFISH-SPECIFIC OVERALL GO AND KEGG ANALYSIS

# Performing GO and KEGG enrichment
run_enrichment <- function(de_raw, prefix, out_dir, plot_dir = NULL) {
  de <- prepare_de(de_raw)
  if (!nrow(de)) return(invisible(NULL))
  gene_mapping <- tryCatch(
    suppressMessages(bitr(unique(de$gene), fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Dr.eg.db)),
    error = function(e) { message(prefix, ": ", conditionMessage(e)); NULL })
  if (is.null(gene_mapping) || !nrow(gene_mapping)) return(invisible(NULL))
  background_ids <- unique(gene_mapping$ENTREZID)
  
  for (direction in c("Upregulated", "Downregulated")) {
    genes <- de$gene[de$regulation == direction]
    ids <- unique(gene_mapping$ENTREZID[gene_mapping$SYMBOL %in% genes])
    if (!length(ids)) next
    
    for (database in c("GO", "KEGG")) {
      result <- tryCatch({
        if (database == "GO") {
          enrichGO(gene = ids, universe = background_ids, OrgDb = org.Dr.eg.db,
                   keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
        } else {
          enrichKEGG(gene = ids, universe = background_ids, organism = "dre",
                     keyType = "ncbi-geneid", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2)
        }
      }, error = function(e) {
        message(prefix, " ", direction, " ", database, ": ", conditionMessage(e))
        NULL
      })
      if (is.null(result)) next
      table_result <- as.data.frame(result)
      
      # Exporting enrichment results
      suffix <- database
      csv_path <- if (is.null(plot_dir)) {
        file.path(out_dir, paste0(database, "_", direction, ".csv"))
      } else {
        file.path(out_dir, database, paste0(prefix, "_", direction, "_", suffix, ".csv"))
      }
      dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(table_result, csv_path, row.names = FALSE)
      if (!nrow(table_result)) next
      
      # Generating enrichment plots
      p <- enrichplot::dotplot(result, showCategory = 15) +
        ggtitle(paste(prefix, direction,
                      if (database == "GO") "GO Enrichment" else "KEGG Pathways")) +
        theme(plot.title = element_text(hjust = 0.5))
      print(p)
      
      # Saving PDF and TIFF
      plot_prefix <- if (is.null(plot_dir)) {
        file.path(out_dir, paste0(database, "_", direction))
      } else {
        file.path(plot_dir, paste0(prefix, "_", direction, "_", suffix))
      }
      dir.create(dirname(plot_prefix), recursive = TRUE, showWarnings = FALSE)
      ggsave(paste0(plot_prefix, ".pdf"), p, width = 11, height = 8)
      ggsave(paste0(plot_prefix, ".tiff"), p, width = 11, height = 8,
             dpi = 300, compression = "lzw")
    }
  }
  invisible(NULL)
}

# Running overall GO and KEGG enrichment
dir.create("Overall_Pathway_Analysis", recursive = TRUE, showWarnings = FALSE)
run_enrichment(overall_DEGs, "Overall", "Overall_Pathway_Analysis")


# 4. CELL TYPE-SPECIFIC GO AND KEGG PATHWAY ANALYSIS

# Creating output directories
for (folder in c("Celltype_Pathway_Analysis/GO", "Celltype_Pathway_Analysis/KEGG",
                 "Celltype_Pathway_Analysis/Plots"))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

# Running cell-type-specific enrichment
for (ct in names(cell_degs)) {
  message("Processing pathway analysis: ", ct)
  run_enrichment(cell_degs[[ct]], ct, "Celltype_Pathway_Analysis",
                 "Celltype_Pathway_Analysis/Plots")
}
message("Cell-type-specific pathway analysis completed.")


# 5. AGGREGATING GO AND KEGG PLOTS

# Creating output folder
dir.create("Integrated_Pathway_Plots", recursive = TRUE, showWarnings = FALSE)

# Reading and combining enrichment results
read_enrichment <- function(database) {
  files <- unlist(lapply(names(cell_degs), function(ct) {
    file.path("Celltype_Pathway_Analysis", database,
              paste0(ct, "_", c("Upregulated", "Downregulated"), "_", database, ".csv"))
  }))
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  results <- lapply(files, function(file) {
    df <- read.csv(file, check.names = FALSE)
    if (!nrow(df)) return(NULL)
    filename <- basename(file)
    df$Celltype <- sub("_(Upregulated|Downregulated)_.*$", "", filename)
    df$Direction <- sub("^.*_(Upregulated|Downregulated)_.*$", "\\1", filename)
    df
  })
  bind_rows(results)
}

# Calculating Jaccard similarity
jaccard_similarity <- function(a, b) {
  denom <- length(union(a, b))
  if (!denom) return(0)
  length(intersect(a, b)) / denom
}

# Creating clustered pathway heatmaps
make_pathway_heatmap <- function(df, database, top_n = 30) {
  if (!nrow(df)) return(NULL)
  df <- df %>% filter(!is.na(p.adjust), !is.na(ID),
                      !is.na(Description), p.adjust < 0.05)
  if (!nrow(df)) return(NULL)
  
  # Selecting top pathways
  top_terms <- df %>% group_by(ID) %>%
    summarise(best_p = min(p.adjust), .groups = "drop") %>%
    slice_min(best_p, n = top_n, with_ties = FALSE) %>% pull(ID)
  df <- df %>% filter(ID %in% top_terms)
  
  # Collecting associated genes
  term_genes <- lapply(top_terms, function(term) {
    unique(unlist(strsplit(paste(df$geneID[df$ID == term], collapse = "/"),
                           "/", fixed = TRUE)))
  })
  
  # Clustering pathways by gene overlap
  similarity_matrix <- outer(seq_along(top_terms), seq_along(top_terms),
                             Vectorize(function(i, j) jaccard_similarity(term_genes[[i]], term_genes[[j]])))
  ordered_terms <- if (length(top_terms) > 1) {
    top_terms[hclust(as.dist(1 - similarity_matrix), method = "average")$order]
  } else top_terms
  
  # Preparing heatmap labels and data
  descriptions <- df %>% distinct(ID, Description)
  labels <- make.unique(descriptions$Description[match(ordered_terms, descriptions$ID)])
  heatmap_df <- df %>%
    mutate(score = -log10(pmax(p.adjust, .Machine$double.xmin))) %>%
    group_by(ID, Celltype, Direction) %>%
    summarise(score = max(score), .groups = "drop") %>%
    mutate(ID = factor(ID, levels = rev(ordered_terms), labels = rev(labels)),
           Celltype = factor(Celltype, levels = celltypes_to_test),
           Direction = factor(Direction, levels = c("Upregulated", "Downregulated")))
  
  # Plotting heatmaps
  p <- ggplot(heatmap_df, aes(Celltype, ID, fill = score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_grid(. ~ Direction) +
    scale_fill_gradient(low = "#DBEAFE", high = "#1D4ED8",
                        name = expression(-log[10]("Adjusted p-value"))) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    theme_minimal(base_size = 12) +
    labs(title = paste(database, "Pathway Enrichment Across Retinal Cell Types"),
         x = "Retinal Cell Type", y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  print(p)
  
  # Saving PDF and TIFF
  plot_prefix <- file.path("Integrated_Pathway_Plots",
                           paste0(database, "_Integrated_Heatmap"))
  plot_height <- max(8, 0.28 * length(top_terms) + 3)
  ggsave(paste0(plot_prefix, ".pdf"), p, width = 24, height = plot_height,
         units = "in", limitsize = FALSE)
  ggsave(paste0(plot_prefix, ".tiff"), p, width = 24, height = plot_height,
         units = "in", limitsize = FALSE, dpi = 300, compression = "lzw")
  p
}

# Generating integrated GO and KEGG heatmaps
p_go_integrated <- make_pathway_heatmap(read_enrichment("GO"), "GO")
p_kegg_integrated <- make_pathway_heatmap(read_enrichment("KEGG"), "KEGG")

sessionInfo()
