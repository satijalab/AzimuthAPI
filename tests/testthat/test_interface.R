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

  annotated <- CloudAzimuth(query, ip = "azimuthapi.satijalab.org")

  expect_s4_class(annotated, "Seurat")
  expect_true(all(c("final_level_labels", "azimuth_label") %in% colnames(annotated@meta.data)))
  expect_identical(unname(annotated@meta.data$azimuth_label), rep(c("T cell", "B cell"), length.out = ncol(query)))
})

test_that("CloudAzimuth builds HTTPS and legacy URLs correctly", {
  expect_identical(
    build_cloud_api_base_url("azimuthapi.satijalab.org"),
    "https://azimuthapi.satijalab.org"
  )

  expect_identical(
    build_cloud_api_base_url("localhost", port = 5000, scheme = "http"),
    "http://localhost:5000"
  )

  expect_identical(
    build_cloud_api_base_url("azimuthapi.satijalab.org", scheme = "http", port = 5000),
    "http://azimuthapi.satijalab.org:5000"
  )

  expect_error(
    build_cloud_api_base_url("azimuthapi.satijalab.org", scheme = "ftp"),
    "`scheme` must be either 'http' or 'https'.",
    fixed = TRUE
  )

})

test_that("CloudAzimuth requires ip without scheme prefix", {
  expect_error(
    CloudAzimuth(object = NULL, ip = "https://example.org"),
    "`ip` should not include 'http://' or 'https://'.",
    fixed = TRUE
  )
})

test_that("CloudAzimuth defaults to official HTTPS endpoint", {
  query <- make_test_object()
  expected <- query
  expected@meta.data$final_level_labels <- rep(c("T cell", "B cell"), length.out = ncol(expected))
  expected@meta.data$azimuth_label <- expected@meta.data$final_level_labels
  api_base_urls <- character()

  testthat::local_mocked_bindings(
    check_api_version = function(url) {
      api_base_urls <<- c(api_base_urls, url)
      invisible(url)
    },
    process_rds_file = function(url, file_path, ...) {
      api_base_urls <<- c(api_base_urls, url)
      saveRDS(expected, file = sub("\\.rds$", "_ANN.rds", file_path))
      invisible(NULL)
    },
    .package = "AzimuthAPI"
  )

  CloudAzimuth(query, ip = "azimuthapi.satijalab.org")
  expect_true(all(api_base_urls[1:2] == "https://azimuthapi.satijalab.org"))
})

test_that("CloudAzimuth preserves alternate-host HTTP defaults", {
  query <- make_test_object()
  expected <- query
  expected@meta.data$final_level_labels <- rep(c("T cell", "B cell"), length.out = ncol(expected))
  expected@meta.data$azimuth_label <- expected@meta.data$final_level_labels
  api_base_urls <- character()

  testthat::local_mocked_bindings(
    check_api_version = function(url) {
      api_base_urls <<- c(api_base_urls, url)
      invisible(url)
    },
    process_rds_file = function(url, file_path, ...) {
      api_base_urls <<- c(api_base_urls, url)
      saveRDS(expected, file = sub("\\.rds$", "_ANN.rds", file_path))
      invisible(NULL)
    },
    .package = "AzimuthAPI"
  )

  CloudAzimuth(query, ip = "localhost")
  expect_true(all(api_base_urls[1:2] == "http://localhost:5000"))
})

test_that("CloudAzimuth preserves legacy official HTTP when port 5000 is explicit", {
  query <- make_test_object()
  expected <- query
  expected@meta.data$final_level_labels <- rep(c("T cell", "B cell"), length.out = ncol(expected))
  expected@meta.data$azimuth_label <- expected@meta.data$final_level_labels
  api_base_urls <- character()

  testthat::local_mocked_bindings(
    check_api_version = function(url) {
      api_base_urls <<- c(api_base_urls, url)
      invisible(url)
    },
    process_rds_file = function(url, file_path, ...) {
      api_base_urls <<- c(api_base_urls, url)
      saveRDS(expected, file = sub("\\.rds$", "_ANN.rds", file_path))
      invisible(NULL)
    },
    .package = "AzimuthAPI"
  )

  CloudAzimuth(query, ip = "azimuthapi.satijalab.org", port = 5000)
  expect_true(all(api_base_urls[1:2] == "http://azimuthapi.satijalab.org:5000"))
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
