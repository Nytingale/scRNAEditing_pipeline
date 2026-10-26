#!/usr/bin/env Rscript
# Install the R packages of the scRNAEditing environment that did not come from
# conda (scRNAEditing.r_packages.tsv), at the exact versions used. Run inside the
# activated environment after creating it from scRNAEditing.linux-64.lock.txt:
#   Rscript envs/install_r_packages.R
# Dependencies are not resolved: every dependency is already in the lock file,
# and resolving would upgrade them.

args    <- commandArgs(trailingOnly = TRUE)
tsv     <- if (length(args)) args[1] else {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  file.path(dirname(normalizePath(f)), "scRNAEditing.r_packages.tsv")
}
pkgs    <- read.delim(tsv, stringsAsFactors = FALSE)
cran    <- "https://cran.r-project.org/src/contrib"
lib     <- .libPaths()[1]
tmp     <- Sys.getenv("TMPDIR", tempdir())

installed_version <- function(p) {
  v <- tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) NA_character_)
  if (is.na(v)) NA_character_ else v
}
same_version <- function(a, b) !is.na(a) && package_version(a) == package_version(b)

install_tarball <- function(url, p) {
  dest <- file.path(tmp, basename(url))
  ok <- tryCatch(utils::download.file(url, dest, quiet = TRUE, mode = "wb") == 0,
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) return(FALSE)
  utils::install.packages(dest, repos = NULL, type = "source", lib = lib,
                          dependencies = FALSE, INSTALL_opts = "--no-test-load")
  unlink(dest)
  TRUE
}

for (i in seq_len(nrow(pkgs))) {
  p <- pkgs$package[i]; v <- pkgs$version[i]; src <- pkgs$source[i]
  if (same_version(installed_version(p), v)) { message(p, " ", v, " already installed"); next }
  message("Installing ", p, " ", v, " (", src, ")")
  if (src == "CRAN") {
    tarball <- sprintf("%s_%s.tar.gz", p, v)
    ok <- install_tarball(file.path(cran, "Archive", p, tarball), p) ||
          install_tarball(file.path(cran, tarball), p)
    if (!ok) stop("Could not download ", p, " ", v, " from CRAN")
  } else if (startsWith(src, "github:")) {
    m <- regmatches(src, regexec("^github:([^@]+)@(.+)$", src))[[1]]
    url <- sprintf("https://codeload.github.com/%s/tar.gz/%s", m[2], m[3])
    dest <- file.path(tmp, paste0(p, "_", m[3], ".tar.gz"))
    utils::download.file(url, dest, quiet = TRUE, mode = "wb")
    exdir <- file.path(tmp, paste0(p, "_src"))
    utils::untar(dest, exdir = exdir)
    utils::install.packages(list.dirs(exdir, recursive = FALSE)[1], repos = NULL,
                            type = "source", lib = lib, dependencies = FALSE)
    unlink(c(dest, exdir), recursive = TRUE)
  } else {
    stop("Unknown source for ", p, ": ", src)
  }
}

got <- vapply(pkgs$package, installed_version, "")
bad <- !mapply(same_version, got, pkgs$version)
if (any(bad)) {
  stop("Version mismatch after install: ",
       paste(sprintf("%s %s (want %s)", pkgs$package[bad], got[bad], pkgs$version[bad]), collapse = ", "))
}
message("All ", nrow(pkgs), " R packages match scRNAEditing.r_packages.tsv")
