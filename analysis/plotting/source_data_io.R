# =============================================================================
# source_data_io.R
# Figure source data: every figure has <stem>_source_data.R, which computes the plotted tables from the
# analysis outputs, and <stem>_plot.R, which reads only those tables and draws the figure.
# Each table is <dir>/<stem>.<part>.tsv plus <stem>.<part>.schema.tsv (column types and factor levels), so the
# plot script gets back exactly the data.table the figure was drawn from (factor order, integer vs double,
# character columns that look numeric).
# =============================================================================
suppressPackageStartupMessages(library(data.table))

# A part is <stem>.<part>.tsv, or <stem>.<part>.tsv.gz when written with gzip = TRUE (large point clouds).
sd_file <- function(dir, stem, part) {
  f <- file.path(dir, sprintf("%s.%s.tsv", stem, part))
  ifelse(!file.exists(f) & file.exists(paste0(f, ".gz")), paste0(f, ".gz"), f)
}
sd_schema_file <- function(dir, stem, part) file.path(dir, sprintf("%s.%s.schema.tsv", stem, part))

# Shortest of 15, 16 or 17 significant digits that reads back as the same double.
exact_double_text <- function(v) {
  out <- sprintf("%.15g", v)
  for (d in 16:17) {
    bad <- which(is.finite(v) & suppressWarnings(as.numeric(out)) != v)
    if (!length(bad)) break
    out[bad] <- sprintf(paste0("%.", d, "g"), v[bad])
  }
  out[is.na(v)] <- NA_character_
  out
}

# Line breaks and tabs in text are stored as the two characters \n and \t (backslash as \\), one record per line.
esc_text <- function(v) gsub("\t", "\\t", gsub("\n", "\\n", gsub("\\", "\\\\", v, fixed = TRUE), fixed = TRUE), fixed = TRUE)
# Text is marked UTF-8 on reading, as R marks "\u2265"-style literals: ggplot2 orders legends by a hash that
# sees the encoding mark.
unesc_text <- function(v) {
  v <- gsub("\\\\", "\001", v, fixed = TRUE)
  v <- gsub("\\t", "\t", gsub("\\n", "\n", v, fixed = TRUE), fixed = TRUE)
  v <- gsub("\001", "\\", v, fixed = TRUE)
  Encoding(v) <- "UTF-8"
  v
}

write_sd <- function(x, dir, stem, part, gzip = FALSE) {
  x <- as.data.table(x)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, sprintf("%s.%s.tsv", stem, part))
  unlink(c(f, paste0(f, ".gz")))
  if (gzip) f <- paste0(f, ".gz")
  schema <- rbindlist(lapply(names(x), function(n) {
    v <- x[[n]]
    type <- if (is.factor(v)) "factor" else if (is.logical(v)) "logical" else if (is.integer(v)) "integer" else
      if (is.numeric(v)) "double" else if (is.character(v)) "character" else
        stop("write_sd: column ", n, " has unsupported class ", class(v)[1])
    data.table(column = n, type = type, level = if (is.factor(v)) levels(v) else NA_character_)
  }))
  for (n in schema[type %in% c("character", "factor"), unique(column)])
    if (any(as.character(x[[n]]) == "NA", na.rm = TRUE)) stop("write_sd: column ", n, " has the text \"NA\", which reads back as missing")
  out <- copy(x)
  for (n in schema[type == "double", column]) set(out, j = n, value = exact_double_text(out[[n]]))
  for (n in schema[type %in% c("character", "factor"), unique(column)])
    set(out, j = n, value = esc_text(as.character(out[[n]])))
  schema[, level := esc_text(level)]
  fwrite(out, f, sep = "\t", na = "NA", quote = "auto", compress = if (gzip) "gzip" else "none")
  fwrite(schema, sd_schema_file(dir, stem, part), sep = "\t", na = "NA", quote = TRUE)
  invisible(f)
}

read_sd <- function(dir, stem, part) {
  f <- sd_file(dir, stem, part)
  if (!file.exists(f)) stop("missing source data ", f)
  schema <- fread(sd_schema_file(dir, stem, part), sep = "\t", na.strings = "NA", colClasses = "character")
  schema[, level := unesc_text(level)]
  cols <- unique(schema$column)
  x <- fread(f, sep = "\t", na.strings = "NA", colClasses = list(character = cols), encoding = "UTF-8")
  stopifnot(identical(names(x), cols))
  for (n in cols) {
    s <- schema[column == n]
    v <- x[[n]]
    if (s$type[1] %in% c("character", "factor")) v <- unesc_text(v)
    v <- switch(s$type[1],
                factor = factor(v, levels = s$level),
                logical = as.logical(v),
                integer = as.integer(v),
                double = as.numeric(v),
                character = v)
    set(x, j = n, value = v)
  }
  x[]
}

# Paths of the tables a plot script read, for the figure's .sources.tsv sidecar.
sd_sources <- function(dir, stem, parts) sd_file(dir, stem, parts)
