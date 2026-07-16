# Note - these tests are designed to test interface and packaging, i.e. they do not interact with the real server to test annotation results

make_test_object <- function() {
  obj <- readRDS(test_path("test_obj.rds"))
  return(Seurat::NormalizeData(obj, verbose = FALSE))
}

test_that("CloudAzimuth returns an annotated Seurat object", {
  query <- make_test_object()
  expected <- query
  expected@meta.data$final_level_labels <- rep(c("T cell", "B cell"), length.out = ncol(expected))
  expected@meta.data$azimuth_label <- expected@meta.data$final_level_labels

  # dummy functions to avoid server calls
  testthat::local_mocked_bindings( 
    check_api_version = function(api_base_url) {
      invisible(api_base_url)
    },
    process_rds_file = function(api_base_url, file_path, ...) { # 
      saveRDS(expected, file = sub("\\.rds$", "_ANN.rds", file_path))
      invisible(NULL)
    },
    .package = "AzimuthAPI"
  )

  annotated <- CloudAzimuth(query, ip = "azimuthapi.satijalab.org", port = 5000)

  expect_s4_class(annotated, "Seurat")
  expect_true(all(c("final_level_labels", "azimuth_label") %in% colnames(annotated@meta.data)))
  expect_identical(unname(annotated@meta.data$azimuth_label), rep(c("T cell", "B cell"), length.out = ncol(query)))
})

test_that("ANNotate returns an annotated Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("reticulate")

  query <- make_test_object()

  # dummy functions to avoid reticulate setup
  fake_import <- function(module, ...) {
    if (identical(module, "panhumanpy.ANNotate")) {
      return(list(
        annotate_core = function(X_query, query_features, cells_meta, annotation_pipeline, eval_batch_size, normalization_override, norm_check_batch_size, 
                                output_mode, refine_labels, map_to_cl, include_cl_id, extract_embeddings, umap_embeddings, n_neighbors, n_components, 
                                metric, min_dist, umap_lr, umap_seed, spread, verbose, init, model_version) {
          cells_meta$final_level_labels <- rep(c("T cell", "B cell"), length.out = nrow(cells_meta))
          list(embeddings_dict = list(), umap_dict = list(), cells_meta = cells_meta)
        }
      ))
    }

    if (identical(module, "scipy.sparse")) {
      return(list(csr_matrix = function(x) x))
    }

    stop(sprintf("Unexpected module request: %s", module))
  }

  testthat::local_mocked_bindings(import = fake_import, r_to_py = function(x) x, .package = "reticulate")

  annotated <- ANNotate(query, extract_embeddings = FALSE, umap_embeddings = FALSE)

  expect_s4_class(annotated, "Seurat")
  expect_true(all(c("final_level_labels", "azimuth_label") %in% colnames(annotated@meta.data)))
  expect_identical(unname(annotated@meta.data$azimuth_label), rep(c("T cell", "B cell"), length.out = ncol(query)))
})
