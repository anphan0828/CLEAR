#!/usr/bin/env Rscript

# Simple runtime and memory benchmark for real-data methods.
# It runs every method script in CLEAR_paper/real_data/methods on 2 TCGA datasets.

get_script_path <- function() {
	cmd_args <- commandArgs(trailingOnly = FALSE)
	file_arg <- grep("^--file=", cmd_args, value = TRUE)
	if (length(file_arg) == 0) {
		return(NA_character_)
	}
	normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

find_project_root <- function(start_dir = getwd()) {
	current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

	repeat {
		has_core <- file.exists(file.path(current, "R", "CLEAR.R"))
		has_real_data <- dir.exists(file.path(current, "CLEAR_paper", "real_data"))
		if (has_core && has_real_data) {
			return(current)
		}

		parent <- dirname(current)
		if (identical(parent, current)) {
			stop("Could not locate the project root containing R/CLEAR.R")
		}
		current <- parent
	}
}

method_args_for <- function(method_file) {
	method_name <- basename(method_file)

	# Shared cutoffs used in the real-data analysis.
	args <- c("-c", "20", "-C", "500")

	# CLEAR variants require iteration/burn-in options.
	if (grepl("^run_CLEAR_", method_name)) {
		args <- c(args, "-i", "1000000", "-b", "500000", "-s", "1")
	}

	args
}

parse_max_rss_kb <- function(stderr_file) {
	if (!file.exists(stderr_file)) {
		return(NA_real_)
	}

	stderr_lines <- readLines(stderr_file, warn = FALSE)
	rss_lines <- grep("^MAX_RSS_KB=", stderr_lines, value = TRUE)
	if (length(rss_lines) == 0) {
		return(NA_real_)
	}

	rss_kb <- suppressWarnings(as.numeric(sub("^MAX_RSS_KB=", "", tail(rss_lines, 1))))
	if (length(rss_kb) == 0 || is.na(rss_kb)) {
		return(NA_real_)
	}

	rss_kb
}

enforce_single_core <- function(n_cores = 1) {
	core_value <- as.character(n_cores)

	thread_env <- c(
		OMP_NUM_THREADS = core_value,
		OMP_THREAD_LIMIT = core_value,
		OPENBLAS_NUM_THREADS = core_value,
		MKL_NUM_THREADS = core_value,
		BLIS_NUM_THREADS = core_value,
		VECLIB_MAXIMUM_THREADS = core_value,
		NUMEXPR_NUM_THREADS = core_value,
		RCPP_PARALLEL_NUM_THREADS = core_value,
		MC_CORES = core_value
	)

	do.call(Sys.setenv, as.list(thread_env))
	options(mc.cores = n_cores)
}

write_single_core_profile <- function(file_path, n_cores = 1) {
	lines <- c(
		paste0("options(mc.cores = ", n_cores, ")"),
		paste0("Sys.setenv(MC_CORES = '", n_cores, "')"),
		"if (requireNamespace('BiocParallel', quietly = TRUE)) {",
		"  BiocParallel::register(BiocParallel::SerialParam())",
		"}"
	)
	writeLines(lines, con = file_path)
}

load_tcga_id_map <- function(dataset_log_file) {
	if (!file.exists(dataset_log_file)) {
		stop(paste("Dataset log file not found:", dataset_log_file))
	}

	if (!requireNamespace("jsonlite", quietly = TRUE)) {
		stop("Package 'jsonlite' is required to read 1_datasets_log.ndjson")
	}

	log_lines <- readLines(dataset_log_file, warn = FALSE)
	log_lines <- log_lines[nzchar(trimws(log_lines))]

	if (length(log_lines) == 0) {
		stop(paste("Dataset log file is empty:", dataset_log_file))
	}

	parsed <- lapply(log_lines, jsonlite::fromJSON)
	is_tcga <- vapply(parsed, function(x) {
		!is.null(x$gene_list) && grepl("/real_data/tcga_data/", x$gene_list)
	}, logical(1))

	parsed_tcga <- parsed[is_tcga]
	if (length(parsed_tcga) == 0) {
		stop("No TCGA entries found in 1_datasets_log.ndjson")
	}

	dataset_names <- vapply(parsed_tcga, function(x) basename(x$gene_list), character(1))
	dataset_ids <- as.integer(vapply(parsed_tcga, function(x) x$ID_dataset, numeric(1)))

	setNames(dataset_ids, dataset_names)
}

run_one_method <- function(method_file,
                        dataset_file,
                        dataset_id,
                        annotation_file,
                        preprocess_script,
                        run_dir,
                        time_exe,
						n_cores) {
	dataset_stem <- tools::file_path_sans_ext(basename(dataset_file))
	method_stem <- tools::file_path_sans_ext(basename(method_file))
	run_tag <- paste(dataset_stem, method_stem, sep = "__")

	result_file <- file.path(run_dir, paste0(run_tag, "_result.tsv"))
	stdout_log <- file.path(run_dir, paste0(run_tag, ".stdout.log"))
	stderr_log <- file.path(run_dir, paste0(run_tag, ".stderr.log"))

	base_args <- c(
		"Rscript", method_file,
		"-g", dataset_file,
		"-a", annotation_file,
		"-o", result_file,
		"-p", preprocess_script,
		method_args_for(method_file)
	)

	start_time <- Sys.time()

	if (file.exists(time_exe)) {
		status <- system2(
			command = time_exe,
			args = c("-f", "MAX_RSS_KB=%M", base_args),
			stdout = stdout_log,
			stderr = stderr_log,
			wait = TRUE
		)
	} else {
		status <- system2(
			command = "Rscript",
			args = base_args[-1],
			stdout = stdout_log,
			stderr = stderr_log,
			wait = TRUE
		)
	}

	elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
	max_rss_kb <- parse_max_rss_kb(stderr_log)

	data.frame(
		dataset_id = dataset_id,
		dataset = basename(dataset_file),
		method = basename(method_file),
		n_cores = n_cores,
		status = ifelse(status == 0, "ok", "failed"),
		exit_code = status,
		elapsed_sec = round(elapsed_sec, 3),
		max_rss_kb = max_rss_kb,
		max_rss_mb = round(max_rss_kb / 1024, 2),
		stringsAsFactors = FALSE
	)
}

main <- function() {
	script_path <- get_script_path()
	start_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
	project_root <- find_project_root(start_dir)
	setwd(project_root)

	methods_dir <- file.path(project_root, "CLEAR_paper", "real_data", "methods")
	tcga_dir <- file.path(project_root, "CLEAR_paper", "real_data", "tcga_data")
	dataset_log_file <- file.path(project_root, "CLEAR_paper", "real_data", "scripts_real_data_analysis", "1_datasets_log.ndjson")
	saved_data_dir <- file.path(project_root, "saved_data")
	bench_root <- file.path(project_root, "CLEAR_paper", "real_data", "results_runtime")
	dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)

	method_files <- sort(list.files(methods_dir, pattern = "^run_.*\\.R$", full.names = TRUE))
	tcga_files <- sort(list.files(tcga_dir, pattern = "\\.tsv$", full.names = TRUE))

	if (length(method_files) == 0) {
		stop("No method scripts found in CLEAR_paper/real_data/methods")
	}
	if (length(tcga_files) < 2) {
		stop("Need at least 2 TCGA datasets in CLEAR_paper/real_data/tcga_data")
	}

	# Simple selection rule: first two TCGA datasets alphabetically.
	selected_tcga <- tcga_files[1:2]
	tcga_id_map <- load_tcga_id_map(dataset_log_file)
	n_cores <- 1

	run_dir <- file.path(bench_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
	dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
	enforce_single_core(n_cores)

	single_core_profile <- file.path(run_dir, "single_core_profile.R")
	write_single_core_profile(single_core_profile, n_cores = n_cores)
	Sys.setenv(R_PROFILE_USER = single_core_profile)

	preprocess_script <- file.path(dirname(script_path), "data_prepare_tcga.R")

	time_exe <- "/usr/bin/time"
	if (!file.exists(time_exe)) {
		warning("/usr/bin/time not found. Runtime will still be recorded but memory may be NA.")
	}
	message("Running benchmark with single-core setting: ", n_cores)

	results <- list()
	idx <- 1

	for (dataset_file in selected_tcga) {
		dataset_name <- basename(dataset_file)
		dataset_id <- as.integer(tcga_id_map[[dataset_name]])
		if (is.na(dataset_id)) {
			warning(paste("Skipping", dataset_name, "because it has no ID mapping in 1_datasets_log.ndjson"))
			next
		}
		annotation_file <- file.path(saved_data_dir, paste0("save_", dataset_id, "_0to100000_annotations.tsv"))

		if (!file.exists(annotation_file)) {
			warning(paste("Skipping", dataset_name, "because annotation file is missing:", annotation_file))
			next
		}

		for (method_file in method_files) {
			message(paste("Running", basename(method_file), "on", dataset_name))
			results[[idx]] <- run_one_method(
				method_file = method_file,
				dataset_file = dataset_file,
				dataset_id = dataset_id,
				annotation_file = annotation_file,
				preprocess_script = preprocess_script,
				run_dir = run_dir,
				time_exe = time_exe,
				n_cores = n_cores
			)
			idx <- idx + 1
		}
	}

	if (length(results) == 0) {
		stop("No method runs were completed.")
	}

	summary_df <- do.call(rbind, results)
	summary_file <- file.path(run_dir, "runtime_memory_summary.tsv")
    # Pretty method names by removing "run_" prefix and ".R" suffix
    summary_df$method <- gsub("^run_|\\.R$", "", summary_df$method)
	write.table(summary_df, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

	message("Saved runtime summary: ", summary_file)
	print(summary_df[, c("dataset", "method", "n_cores", "status", "elapsed_sec", "max_rss_mb")], row.names = FALSE)
}

main()
