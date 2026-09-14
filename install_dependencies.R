#!/usr/bin/env Rscript

# install_dependencies.R
# Run this script before installing multiBAT.combo:
#   Rscript install_dependencies.R
#
# Place this file in the package root directory, next to DESCRIPTION.

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Whether to install packages listed only under Suggests.
# Set to TRUE if you want to run examples/tests that require suggested packages.
install_suggests <- FALSE

# Locate this script and the package DESCRIPTION file.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  package_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
} else {
  package_dir <- getwd()
}

description_file <- file.path(package_dir, "DESCRIPTION")

if (!file.exists(description_file)) {
  stop(
    "DESCRIPTION was not found in:\n  ", package_dir,
    "\n\nPlace install_dependencies.R in the package root directory ",
    "(the same directory as DESCRIPTION) and run:\n",
    "  Rscript install_dependencies.R"
  )
}

cat("Package directory:\n  ", package_dir, "\n\n", sep = "")

desc <- read.dcf(description_file)

fields <- c("Depends", "Imports", "LinkingTo")
if (isTRUE(install_suggests)) {
  fields <- c(fields, "Suggests")
}

# Collect dependency specifications from DESCRIPTION.
dep_text <- unlist(desc[1, intersect(fields, colnames(desc)), drop = TRUE])
dep_text <- dep_text[!is.na(dep_text) & nzchar(dep_text)]

if (length(dep_text) == 0) {
  cat("No package dependencies were found in DESCRIPTION.\n")
  quit(status = 0)
}

dep_entries <- trimws(unlist(strsplit(paste(dep_text, collapse = ","), ",")))
dep_entries <- dep_entries[nzchar(dep_entries)]

parse_dependency <- function(x) {
  pkg <- trimws(sub("\\s*\\(.*\\)\\s*$", "", x))
  constraint <- NA_character_

  if (grepl("\\(", x)) {
    constraint <- sub("^.*\\((.*)\\).*$", "\\1", x)
    constraint <- trimws(constraint)
  }

  data.frame(
    package = pkg,
    constraint = constraint,
    stringsAsFactors = FALSE
  )
}

deps <- do.call(rbind, lapply(dep_entries, parse_dependency))
deps <- deps[deps$package != "R" & nzchar(deps$package), , drop = FALSE]
deps <- deps[!duplicated(deps$package), , drop = FALSE]

if (nrow(deps) == 0) {
  cat("No external R package dependencies were found in DESCRIPTION.\n")
  quit(status = 0)
}

version_ok <- function(installed_version, constraint) {
  if (is.na(constraint) || !nzchar(constraint)) {
    return(TRUE)
  }

  m <- regexec("^\\s*(>=|<=|==|=|>|<)\\s*([^[:space:]]+)\\s*$", constraint)
  parts <- regmatches(constraint, m)[[1]]

  # If a constraint cannot be parsed, leave version validation to R during
  # package installation rather than failing here.
  if (length(parts) != 3) {
    warning("Could not parse version constraint: ", constraint)
    return(TRUE)
  }

  op <- parts[2]
  required_version <- parts[3]
  cmp <- compareVersion(as.character(installed_version), required_version)

  switch(
    op,
    ">=" = cmp >= 0,
    "<=" = cmp <= 0,
    "==" = cmp == 0,
    "="  = cmp == 0,
    ">"  = cmp > 0,
    "<"  = cmp < 0,
    TRUE
  )
}

installed <- installed.packages()
to_install <- character(0)

cat("Checking R package dependencies...\n\n")

for (i in seq_len(nrow(deps))) {
  pkg <- deps$package[i]
  constraint <- deps$constraint[i]

  if (!(pkg %in% rownames(installed))) {
    cat(sprintf("[MISSING] %s", pkg))
    if (!is.na(constraint)) cat(sprintf(" (%s)", constraint))
    cat("\n")
    to_install <- c(to_install, pkg)
    next
  }

  current_version <- packageVersion(pkg)

  if (!version_ok(current_version, constraint)) {
    cat(sprintf(
      "[UPDATE ] %s %s installed; required %s\n",
      pkg, as.character(current_version), constraint
    ))
    to_install <- c(to_install, pkg)
  } else {
    cat(sprintf("[OK     ] %s %s\n", pkg, as.character(current_version)))
  }
}

to_install <- unique(to_install)

if (length(to_install) == 0) {
  cat("\nAll required R package dependencies are already installed.\n")
  quit(status = 0)
}

cat(
  "\nInstalling missing or outdated packages from CRAN:\n  ",
  paste(to_install, collapse = ", "),
  "\n\n",
  sep = ""
)

install.packages(
  to_install,
  dependencies = c("Depends", "Imports", "LinkingTo")
)

# Verify installation.
failed <- character(0)

for (i in seq_len(nrow(deps))) {
  pkg <- deps$package[i]
  constraint <- deps$constraint[i]

  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
    next
  }

  if (!version_ok(packageVersion(pkg), constraint)) {
    failed <- c(failed, pkg)
  }
}

failed <- unique(failed)

if (length(failed) > 0) {
  stop(
    "\nThe following dependencies could not be installed or do not meet ",
    "the required version:\n  ",
    paste(failed, collapse = ", "),
    "\n\nCheck the messages above. If any dependency is hosted outside CRAN ",
    "(for example, Bioconductor or GitHub), install it separately and rerun ",
    "this script."
  )
}

cat("\nAll required dependencies are installed successfully.\n")
cat("You can now install the package, for example:\n")
cat("  R CMD INSTALL .\n")
