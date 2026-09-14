#!/usr/bin/env Rscript

# install_dependencies.R
#
# Run before installing multiBAT.combo:
#
#   Rscript multiBAT.combo-main/install_dependencies.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_suggests <- FALSE

# ------------------------------------------------------------
# Locate DESCRIPTION
# ------------------------------------------------------------

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
    "\n\nPlease place install_dependencies.R in the package root directory ",
    "(the same directory as DESCRIPTION)."
  )
}

cat("Package directory:\n  ", package_dir, "\n\n", sep = "")

desc <- read.dcf(description_file)

# ------------------------------------------------------------
# Check R version
# ------------------------------------------------------------

if ("Depends" %in% colnames(desc)) {

  depends_text <- desc[1, "Depends"]

  if (!is.na(depends_text) && grepl("R\\s*\\(", depends_text)) {

    r_entry <- trimws(
      grep(
        "^R\\s*\\(",
        trimws(unlist(strsplit(depends_text, ","))),
        value = TRUE
      )
    )

    if (length(r_entry) > 0) {

      constraint <- sub("^R\\s*\\((.*)\\)$", "\\1", r_entry)
      constraint <- trimws(constraint)

      m <- regexec(
        "^\\s*(>=|<=|==|=|>|<)\\s*([^[:space:]]+)\\s*$",
        constraint
      )

      parts <- regmatches(constraint, m)[[1]]

      if (length(parts) == 3) {

        op <- parts[2]
        required_version <- parts[3]

        cmp <- compareVersion(
          as.character(getRversion()),
          required_version
        )

        r_ok <- switch(
          op,
          ">=" = cmp >= 0,
          "<=" = cmp <= 0,
          "==" = cmp == 0,
          "="  = cmp == 0,
          ">"  = cmp > 0,
          "<"  = cmp < 0
        )

        if (!r_ok) {
          stop(
            "R version requirement not satisfied.\n",
            "Installed: ", getRversion(), "\n",
            "Required: ", constraint
          )
        }

        cat(
          "[OK     ] R ",
          getRversion(),
          " (required ",
          constraint,
          ")\n\n",
          sep = ""
        )
      }
    }
  }
}

# ------------------------------------------------------------
# Read package dependencies
# ------------------------------------------------------------

fields <- c("Depends", "Imports", "LinkingTo")

if (isTRUE(install_suggests)) {
  fields <- c(fields, "Suggests")
}

fields <- intersect(fields, colnames(desc))

dep_text <- unlist(desc[1, fields, drop = TRUE])
dep_text <- dep_text[!is.na(dep_text) & nzchar(dep_text)]

if (length(dep_text) == 0) {
  cat("No package dependencies were found in DESCRIPTION.\n")
  quit(status = 0)
}

dep_entries <- trimws(
  unlist(strsplit(paste(dep_text, collapse = ","), ","))
)

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

deps <- do.call(
  rbind,
  lapply(dep_entries, parse_dependency)
)

# Remove R itself
deps <- deps[
  deps$package != "R" & nzchar(deps$package),
  ,
  drop = FALSE
]

deps <- deps[
  !duplicated(deps$package),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Version checking
# ------------------------------------------------------------

version_ok <- function(installed_version, constraint) {

  if (is.na(constraint) || !nzchar(constraint)) {
    return(TRUE)
  }

  m <- regexec(
    "^\\s*(>=|<=|==|=|>|<)\\s*([^[:space:]]+)\\s*$",
    constraint
  )

  parts <- regmatches(constraint, m)[[1]]

  if (length(parts) != 3) {
    warning("Could not parse version constraint: ", constraint)
    return(TRUE)
  }

  op <- parts[2]
  required_version <- parts[3]

  cmp <- compareVersion(
    as.character(installed_version),
    required_version
  )

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

# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

installed <- installed.packages()
to_install <- character(0)

cat("Checking R package dependencies...\n\n")

for (i in seq_len(nrow(deps))) {

  pkg <- deps$package[i]
  constraint <- deps$constraint[i]

  if (!(pkg %in% rownames(installed))) {

    cat(sprintf("[MISSING] %s", pkg))

    if (!is.na(constraint)) {
      cat(sprintf(" (%s)", constraint))
    }

    cat("\n")

    to_install <- c(to_install, pkg)
    next
  }

  current_version <- packageVersion(pkg)

  if (!version_ok(current_version, constraint)) {

    cat(
      sprintf(
        "[UPDATE ] %s %s installed; required %s\n",
        pkg,
        as.character(current_version),
        constraint
      )
    )

    to_install <- c(to_install, pkg)

  } else {

    cat(
      sprintf(
        "[OK     ] %s %s\n",
        pkg,
        as.character(current_version)
      )
    )
  }
}

to_install <- unique(to_install)

# ------------------------------------------------------------
# Install missing dependencies
# ------------------------------------------------------------

if (length(to_install) == 0) {

  cat(
    "\nAll required R package dependencies are already installed.\n"
  )

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
  dependencies = c(
    "Depends",
    "Imports",
    "LinkingTo"
  )
)

# ------------------------------------------------------------
# Verify installation
# ------------------------------------------------------------

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
    "\nThe following dependencies could not be installed ",
    "or do not meet the required version:\n  ",
    paste(failed, collapse = ", "),
    "\n\nPlease check the installation messages above."
  )
}

cat(
  "\nAll required dependencies were installed successfully.\n"
)

cat(
  "\nYou can now install multiBAT.combo.\n"
)
