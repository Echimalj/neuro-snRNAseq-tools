#' Save object checkpoint
#'
#' @param object R object to save.
#' @param file Output RDS file.
#'
#' @return Invisibly returns file path.
#' @export
save_checkpoint <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  saveRDS(object, file = file)

  message("Saved checkpoint: ", file)
  invisible(file)
}

#' Load object checkpoint
#'
#' @param file RDS file to load.
#'
#' @return Loaded R object.
#' @export
load_checkpoint <- function(file) {
  if (!file.exists(file)) {
    stop("Checkpoint file does not exist: ", file, call. = FALSE)
  }

  message("Loading checkpoint: ", file)
  readRDS(file)
}

#' Save checkpoint and optionally reload it
#'
#' Useful after memory-heavy steps because reloading the object can release
#' memory from intermediate objects.
#'
#' @param object R object to save.
#' @param file Output RDS file.
#' @param reload Logical. If TRUE, reloads the object after saving.
#'
#' @return Original or reloaded object.
#' @export
checkpoint <- function(object, file, reload = TRUE) {
  save_checkpoint(object, file)

  rm(object)
  gc()

  if (reload) {
    return(load_checkpoint(file))
  }

  invisible(NULL)
}
