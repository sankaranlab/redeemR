#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    input_rds = NULL,
    output_html = NULL,
    sample_name = NULL,
    thr = NULL,
    raw_name = 'RawGenotypes.Sensitive.StrandBalance',
    top_n_variants = 20L,
    consensus_n = 25L,
    output_metrics_tsv = NULL,
    output_report_rds = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, '--')) {
      stop('Unexpected argument: ', key, call. = FALSE)
    }
    val <- if (i < length(args)) args[[i + 1L]] else NULL
    switch(
      key,
      '--input-rds' = out$input_rds <- val,
      '--output-html' = out$output_html <- val,
      '--sample-name' = out$sample_name <- val,
      '--thr' = out$thr <- val,
      '--raw-name' = out$raw_name <- val,
      '--top-n-variants' = out$top_n_variants <- as.integer(val),
      '--consensus-n' = out$consensus_n <- as.integer(val),
      '--output-metrics-tsv' = out$output_metrics_tsv <- val,
      '--output-report-rds' = out$output_report_rds <- val,
      stop('Unknown argument: ', key, call. = FALSE)
    )
    i <- i + 2L
  }

  if (is.null(out$input_rds) || is.null(out$output_html)) {
    stop(
      paste(
        'Usage:',
        'render_filter2_qc.R --input-rds <path> --output-html <path>',
        '[--sample-name <name>] [--thr <thr>]',
        '[--raw-name <raw file>] [--top-n-variants <n>]',
        '[--consensus-n <n>] [--output-metrics-tsv <path>]',
        '[--output-report-rds <path>]'
      ),
      call. = FALSE
    )
  }

  out
}

first_existing <- function(paths, dir = FALSE) {
  for (p in paths) {
    if (is.na(p) || is.null(p) || !nzchar(p)) next
    if (dir) {
      if (dir.exists(p)) return(p)
    } else {
      if (file.exists(p)) return(p)
    }
  }
  NULL
}

move_if_present <- function(src_candidates, dest, dir = FALSE) {
  src <- first_existing(src_candidates, dir = dir)
  if (is.null(src)) return(FALSE)
  if (dir) {
    if (dir.exists(dest) || file.exists(dest)) unlink(dest, recursive = TRUE, force = TRUE)
  } else {
    if (file.exists(dest)) unlink(dest, force = TRUE)
  }
  ok <- file.rename(src, dest)
  isTRUE(ok)
}

script_path <- normalizePath(
  sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]),
  mustWork = TRUE
)
repo_root <- dirname(dirname(dirname(script_path)))
template_path <- file.path(repo_root, 'inst', 'quarto', 'filter2_qc_report.qmd')
asset_dir_name <- 'filter2_qc_report_files'

if (!file.exists(template_path)) {
  stop('Cannot find Quarto template at ', template_path, call. = FALSE)
}

opts <- parse_args(args)
process_cwd <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- normalizePath(dirname(opts$output_html), mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
final_html <- normalizePath(opts$output_html, mustWork = FALSE)
final_asset_dir <- file.path(output_dir, asset_dir_name)
staging_dir <- file.path(output_dir, paste0('.render_', basename(opts$output_html)))
if (dir.exists(staging_dir) || file.exists(staging_dir)) {
  unlink(staging_dir, recursive = TRUE, force = TRUE)
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
staging_html <- file.path(staging_dir, basename(opts$output_html))
staging_asset_dir <- file.path(staging_dir, asset_dir_name)
local_template_path <- file.path(staging_dir, basename(template_path))
invisible(file.copy(template_path, local_template_path, overwrite = TRUE))

render_with_quarto_pkg <- function() {
  if (!requireNamespace('quarto', quietly = TRUE)) {
    return(FALSE)
  }
  quarto::quarto_render(
    input = local_template_path,
    output_file = basename(opts$output_html),
    output_dir = staging_dir,
    execute_params = list(
      input_rds = opts$input_rds,
      sample_name = opts$sample_name,
      thr = opts$thr,
      raw_name = opts$raw_name,
      top_n_variants = opts$top_n_variants,
      consensus_n = opts$consensus_n,
      output_metrics_tsv = if (is.null(opts$output_metrics_tsv)) NULL else file.path(staging_dir, basename(opts$output_metrics_tsv)),
      output_report_rds = if (is.null(opts$output_report_rds)) NULL else file.path(staging_dir, basename(opts$output_report_rds))
    ),
    quiet = FALSE
  )
  TRUE
}

render_with_cli <- function() {
  cli <- Sys.which('quarto')
  if (!nzchar(cli)) {
    stop('Quarto CLI not found and quarto R package is unavailable.', call. = FALSE)
  }
  qargs <- c(
    'render', local_template_path,
    '--to', 'html',
    '--output', basename(staging_html),
    '--output-dir', staging_dir,
    '-P', paste0('input_rds:', normalizePath(opts$input_rds)),
    '-P', paste0('sample_name:', ifelse(is.null(opts$sample_name), '', opts$sample_name)),
    '-P', paste0('thr:', ifelse(is.null(opts$thr), '', opts$thr)),
    '-P', paste0('raw_name:', opts$raw_name),
    '-P', paste0('top_n_variants:', opts$top_n_variants),
    '-P', paste0('consensus_n:', opts$consensus_n),
    '-P', paste0('output_metrics_tsv:', ifelse(is.null(opts$output_metrics_tsv), '', file.path(staging_dir, basename(opts$output_metrics_tsv)))),
    '-P', paste0('output_report_rds:', ifelse(is.null(opts$output_report_rds), '', file.path(staging_dir, basename(opts$output_report_rds))))
  )
  status <- system2(cli, args = qargs)
  if (!identical(status, 0L)) {
    stop('Quarto CLI render failed with exit status ', status, call. = FALSE)
  }
}

if (!render_with_quarto_pkg()) {
  render_with_cli()
}

rendered_html_candidates <- c(
  staging_html,
  file.path(process_cwd, basename(opts$output_html)),
  file.path(output_dir, basename(opts$output_html))
)
rendered_asset_candidates <- c(
  staging_asset_dir,
  file.path(process_cwd, asset_dir_name),
  file.path(output_dir, asset_dir_name)
)

if (!move_if_present(rendered_html_candidates, final_html, dir = FALSE)) {
  stop('Failed to locate and move rendered HTML into final location.', call. = FALSE)
}
move_if_present(rendered_asset_candidates, final_asset_dir, dir = TRUE)

if (!is.null(opts$output_metrics_tsv) && nzchar(opts$output_metrics_tsv)) {
  metrics_candidates <- c(
    file.path(staging_dir, basename(opts$output_metrics_tsv)),
    file.path(process_cwd, basename(opts$output_metrics_tsv)),
    file.path(output_dir, basename(opts$output_metrics_tsv))
  )
  move_if_present(metrics_candidates, opts$output_metrics_tsv, dir = FALSE)
}

if (!is.null(opts$output_report_rds) && nzchar(opts$output_report_rds)) {
  report_candidates <- c(
    file.path(staging_dir, basename(opts$output_report_rds)),
    file.path(process_cwd, basename(opts$output_report_rds)),
    file.path(output_dir, basename(opts$output_report_rds))
  )
  move_if_present(report_candidates, opts$output_report_rds, dir = FALSE)
}

qc_png_name <- paste0(
  ifelse(is.null(opts$sample_name) || !nzchar(opts$sample_name), sub('[.].*$', '', basename(opts$input_rds)), opts$sample_name),
  '.redeemV_qc.png'
)
qc_png_dest <- file.path(output_dir, qc_png_name)
qc_png_candidates <- c(
  file.path(staging_dir, qc_png_name),
  file.path(process_cwd, qc_png_name),
  file.path(output_dir, qc_png_name)
)
move_if_present(qc_png_candidates, qc_png_dest, dir = FALSE)

make_self_contained <- function() {
  cli <- Sys.which('quarto')
  if (!nzchar(cli)) return(FALSE)
  tmp_html <- paste0(basename(final_html), '.selfcontained')
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(output_dir)
  status <- system2(
    cli,
    args = c(
      'pandoc',
      basename(final_html),
      '--from', 'html',
      '--to', 'html',
      '--embed-resources',
      '--standalone',
      '-o', tmp_html
    )
  )
  if (!identical(status, 0L) || !file.exists(tmp_html)) {
    return(FALSE)
  }
  unlink(basename(final_html), force = TRUE)
  ok <- file.rename(tmp_html, basename(final_html))
  if (!isTRUE(ok)) {
    return(FALSE)
  }
  if (dir.exists(final_asset_dir) || file.exists(final_asset_dir)) {
    unlink(final_asset_dir, recursive = TRUE, force = TRUE)
  }
  if (file.exists(qc_png_dest)) {
    unlink(qc_png_dest, force = TRUE)
  }
  TRUE
}

invisible(make_self_contained())

cleanup_candidates <- c(
  file.path(process_cwd, basename(opts$output_html)),
  file.path(process_cwd, asset_dir_name),
  file.path(process_cwd, qc_png_name)
)
for (p in cleanup_candidates) {
  if (dir.exists(p)) unlink(p, recursive = TRUE, force = TRUE)
  if (file.exists(p)) unlink(p, force = TRUE)
}

unlink(staging_dir, recursive = TRUE, force = TRUE)
