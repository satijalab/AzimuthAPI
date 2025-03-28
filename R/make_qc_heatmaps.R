library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

make_QC_heatmap <- function(
    seurat_obj, 
    group.by=NULL, 
    cells.order=NULL, 
    n_downsample = 200, 
    save_folder_path = NULL, 
    logfc_cutoff = log(2), 
    n_markers = 10, 
    text.size=5, 
    text.angle=90, 
    min.size = NULL, 
    max.size = NULL, 
    min.pct=0.1, 
    reorder = TRUE, 
    switch_id=NULL,
    identity = "" 
) {
  if (!is.null(group.by)) Idents(seurat_obj) <- group.by
  if (!is.null(min.size)) {
    seurat_obj <- subset(seurat_obj,idents = names(which(table(Idents(seurat_obj))>=min.size)))
  }
  if (!is.null(max.size)) {
    seurat_obj <- subset(seurat_obj,idents = names(which(table(Idents(seurat_obj))<=max.size)))
  }
  if (!is.null(n_downsample)) seurat_obj <- subset(seurat_obj, downsample=n_downsample)
  mark_all <- FindAllMarkers(seurat_obj,only.pos = TRUE,min.pct = min.pct)
  mark_all %>% group_by(cluster) %>%
    dplyr::filter(avg_log2FC > logfc_cutoff) %>%
    slice_head(n = n_markers) %>%
    ungroup() -> top_markers
  
  if (!is.null(switch_id)) {
    Idents(seurat_obj) <- switch_id
  }
  
  seurat_obj <- ScaleData(seurat_obj,features = top_markers$gene)

  if (reorder) {
    # Compute average profiles for each cell type
    average_profiles <- LayerData(seurat_obj,layer = 'scale.data') %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Cell") %>%
      left_join(data.frame(Cell = Cells(seurat_obj), CellType = Idents(seurat_obj)), by = "Cell") %>%
      group_by(CellType) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% # Use `where(is.numeric)` to avoid non-numeric columns
      column_to_rownames("CellType")
    

    dist_matrix <- dist(as.matrix(cor(t(average_profiles))))
    hclust_results <- hclust(dist_matrix)
    
    # Extract the ordered cell types from the dendrogram
    ordered_cell_types <- hclust_results$labels[hclust_results$order]
    
    # Step 4: Reorder the Idents of the Seurat object based on the dendrogram
    seurat_obj <- SetIdent(seurat_obj, value = factor(Idents(seurat_obj), levels = ordered_cell_types))
    
    # Step 5: Reorder the top_markers dataframe based on ordered_cell_types
    top_markers <- top_markers %>%
      mutate(cluster = factor(cluster, levels = ordered_cell_types)) %>%
      arrange(cluster)
  }
  
  cells.plot <- names(which(!is.na(Idents(seurat_obj))))
  if (!is.null(cells.order)) cells.plot <- intersect(cells.order,cells.plot)
  plot_heatmap <- DoHeatmap(seurat_obj,features = top_markers$gene,cells = cells.plot, size = 3.5, angle = text.angle)+theme(
    axis.text.y = element_text(size = text.size)) + NoLegend()
  if (!is.null(save_folder_path))  ggsave(paste0(save_folder_path, identity, "_heatmap.png"))
  return(plot_heatmap)
}

make_azimuth_QC_heatmaps <- function(
    object, 
    final_name = 'azimuth_fine', 
    level1_name = 'azimuth_broad', 
    full_name = 'full_hierarchical_labels', 
    min.final.group = 10, 
    max.ids.per.plot = 10, 
    ...
) {
  # Ensure the required columns exist
  required_columns <- c(level1_name,full_name,final_name)
  if (!all(required_columns %in% colnames(object@meta.data))) { #SKYLAR: Proposed fix
    stop("The Seurat object is missing required metadata.")
  }
  
  # Extract metadata
  metadata <- object@meta.data
  
  # filter non-concordant
  # SKYLAR: Proposed fix due to ANNotate metadata column name change
  metadata <- subset(metadata, full_consistent_hierarchy==TRUE)

  level2_name = 'level_two_labels'
  metadata[, level2_name] <- sapply(strsplit(metadata[, full_name], "\\|"), function(x) ifelse(length(x) > 1, x[2], ""))
  
  abundance_filter_pass <- names(which(table(metadata[,final_name]) > min.final.group))
  metadata <- metadata[which(metadata[,final_name] %in% abundance_filter_pass),]
  
  # Adjust first-level grouping: Separate "Immune cell" into "Immune_Lymphoid" and "Immune_Myeloid"
  metadata$adjusted_level_1 <- ifelse(
    metadata[,level1_name] == "Immune cell" & metadata[,level2_name] %in% c("Myeloid cell","Lymphoid cell"),
    paste0("Immune_", metadata[,level2_name]),
    as.character(metadata[,level1_name]) #SKYLAR: Proposed fix
  )
  
  # anything that doesn't fit (i.e. erythroid) goes into myeloid category
  metadata[metadata$adjusted_level_1=="Immune cell","adjusted_level_1"] <- 'Immune_Myeloid cell'
  
  # Split by adjusted first level
  split_list <- split(metadata, metadata$adjusted_level_1)
  
  result_list <- list()
  
  for (level1 in names(split_list)) {
    subset_df <- split_list[[level1]]
    
    # Split by second level
    second_level_groups <- split(subset_df, subset_df[,final_name])
    
    # Remove groups that do not meet the min.final.group threshold
    second_level_groups <- second_level_groups[sapply(second_level_groups, nrow) >= min.final.group]
    
    # Check if more than 10 groups and split further
    second_level_names <- names(sort(unlist(lapply(second_level_groups,nrow)),decreasing = TRUE))
    num_levels <- length(second_level_names)
    
    if (is.null(max.ids.per.plot) || num_levels <= max.ids.per.plot) {
        if (num_levels > 0) {
            group_name <- paste0(level1, "_1") # If only one group, assign "_1" suffix
            result_list[[group_name]] <- subset_df
        }
    } else {
      num_splits <- ceiling(num_levels / max.ids.per.plot)
      split_indices <- rep(1:num_splits, each = max.ids.per.plot, length.out = num_levels)
      
      for (split_idx in unique(split_indices)) {
        selected_labels <- second_level_names[split_indices == split_idx]
        group_name <- paste0(level1, "_", split_idx) # Naming convention: "level1_1", "level1_2", etc.
        result_list[[group_name]] <- subset_df[subset_df[,final_name] %in% selected_labels, ]
      }
    } 
  }

   if(!is.null(cells.order) && length(cells.order) == 1){
       cells.order <- rownames(seurat_obj[[]])[order(seurat_obj[[]][[cells.order]])]
   } 

  plot_list <- list()
  for (level1 in names(result_list)) {
    lobj <- subset(object,cells = rownames(result_list[[level1]]))
    Idents(lobj) <- final_name
    tryCatch({
      plot_list[[level1]] <- make_QC_heatmap(lobj, min.size = min.final.group, identity = as.character(level1), ...)
    }, error = function(e) {
      message(paste("Error in processing", level1, ":", e$message))
    })
  }
  return(plot_list)
}
