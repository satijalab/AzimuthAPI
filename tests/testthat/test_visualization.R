
make_viz_test_object <- function(cells = NULL) {
  query <- readRDS(test_path("test_obj.rds"))
  query <- Seurat::NormalizeData(query, verbose = FALSE)
  if (!is.null(cells)) {
    query <- subset(query, cells = colnames(query)[seq_len(cells)])
  }
  return(query)
}

set_azimuth_metadata <- function(query, fine, broad, full) {
  query@meta.data$azimuth_fine <- fine
  query@meta.data$azimuth_broad <- broad
  query@meta.data$full_hierarchical_labels <- full
  query@meta.data$full_consistent_hierarchy <- TRUE
  return(query)
}

test_that("make_QC_heatmap rejects single-label inputs before marker detection", {
  query <- make_viz_test_object()
  query@meta.data$final_level_labels <- "Unassigned"

  expect_error(
    make_QC_heatmap(query, group.by = "final_level_labels"),
    "at least two identity groups",
    fixed = TRUE
  )
})

test_that("make_azimuth_QC_heatmaps skips single-label final annotation groups", {
  query <- make_viz_test_object()
  query <- set_azimuth_metadata(
    query,
    fine = "Unassigned",
    broad = "Unassigned",
    full = "Unassigned"
  )

  expect_warning(
    plots <- make_azimuth_QC_heatmaps(query, min.final.group = 10),
    "Skipping Unassigned_1",
    fixed = TRUE
  )
  expect_false("Unassigned_1" %in% names(plots))
})

test_that("make_azimuth_QC_heatmaps returns valid plots even when one group errors", {
  query <- make_viz_test_object(cells = 70)
  query <- set_azimuth_metadata(
    query,
    fine = c(
      rep("Unassigned", 10),
      rep("T cell", 10),
      rep("B cell", 10),
      rep("Fibroblast", 25),
      rep("Rare stromal", 15)
    ),
    broad = c(
      rep("Immune cell", 30),
      rep("Stromal cell", 40)
    ),
    full = c(
      rep("Immune cell|Lymphoid cell", 30),
      rep("Stromal cell", 40)
    )
  )

  # after the filtering, the stromal cell group has only "Rare stromal" left, but this is a single-label group, so it will error out in make_QC_heatmap
  # the valid plot for "Immune cell" should still be returned, and the error for "Stromal cell" should be caught and reported as a warning
  expect_warning(
    plots <- make_azimuth_QC_heatmaps(
      query,
      min.final.group = 9,
      n_markers = 1,
      max.size = 20
    ),
    "Error in processing Stromal cell_1", 
    fixed = TRUE
  )

  expect_named(plots, "Immune_Lymphoid cell_1")
  expect_s3_class(plots[["Immune_Lymphoid cell_1"]], "ggplot")
  expect_false("Stromal cell_1" %in% names(plots))
})
