# Utilities 

find_data_file <- function(base_dir, filename) {
  # Recursive search for a file in base_dir
  files <- list.files(base_dir, pattern = filename, recursive = TRUE, full.names = TRUE)
  if (length(files) > 0) {
    return(files[1])  # return first match
  } else {
    return(NULL)
  }
}
