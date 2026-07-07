
check_file_exist <- function(basename) {
  if(missing(basename)) return
  if(basename == "") return
  hbd.file  <- paste0(path.expand(basename), ".hbd")
  flod.file <- paste0(path.expand(basename), ".flod")
  if(file.exists(hbd.file) || file.exists(flod.file)) {
    stop("File already exist\n")
  }
}

