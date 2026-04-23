#' Check that required packages are installed
#'
#' @param packages Character vector of package names.
#'
#' @return Invisibly returns TRUE if all packages are installed.
#' @export
check_required_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this workflow.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
