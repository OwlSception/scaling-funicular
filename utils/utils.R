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

def update_requirements(project_directory):
    """Updates the requirements.txt file with the current packages."""
    try:
        req_path = os.path.join(project_directory, "docs", "requirements.txt")
        with open(req_path, "w") as req_file:
            subprocess.run(["pip", "freeze"], stdout=req_file)
        print(f'Updated requirements in {req_path}')
    except Exception as e:
        print(f'Error updating requirements: {e}')