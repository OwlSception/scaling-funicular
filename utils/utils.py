# Utilities and Functions 

from pathlib import Path

def find_data_file(base_dir, filename):
    """
    Search for a file within a directory and its subdirectories.

    Args:
        base_dir (str or Path): The root folder to start searching from.
        filename (str): The exact filename to find (e.g., "data.csv").

    Returns:
        Path or None: Full path to the first match, or None if not found.
    """
    base_dir = Path(base_dir)
    for file_path in base_dir.rglob(filename):  # recursive search
        return file_path
    return None

# Example usage:
file_path = find_data_file("raw_data", "data.csv")
if file_path:
    print("Found file at:", file_path)
else:
    print("File not found")

def update_requirements(project_directory):
    """Updates the requirements.txt file with the current packages."""
    try:
        req_path = os.path.join(project_directory, "docs", "requirements.txt")
        with open(req_path, "w") as req_file:
            subprocess.run(["pip", "freeze"], stdout=req_file)
        print(f'Updated requirements in {req_path}')
    except Exception as e:
        print(f'Error updating requirements: {e}')