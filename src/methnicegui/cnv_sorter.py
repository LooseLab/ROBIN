import os
import shutil
import tempfile
from natsort import natsorted


class SymbolicLinkCreator:
    """
    Class to create a subset of symbolic links to files in a given folder.
    This will be used to help handle multiple bam files to recreate copy number plots over time.
    """

    def __init__(self, bam_path):
        self.bam_path = bam_path
        self.temp_folder = None
        pass

    def _create_temporary_folder(self):
        """
        Create a temporary folder to store the symbolic links.
        """
        if not self.temp_folder:
            self.temp_folder = tempfile.mkdtemp()

    def iterate_bam_files(self):
        """
        Iterate over the bam files in the given folder.
        """

        counter = 0
        bam_list = []
        for file in natsorted(os.listdir(os.path.abspath(self.bam_path))):
            if file.endswith(".bam"):
                counter += 1
                print(os.path.join(self.bam_path, file))
                bam_list.append(os.path.join(self.bam_path, file))
                os.mkdir(os.path.join(self.temp_folder, f"{counter}"))
                for file_path in bam_list:
                    file_name = os.path.basename(file_path)
                    link_path = os.path.join(self.temp_folder, f"{counter}", file_name)
                    os.symlink(file_path, link_path)
                    print(f"Created symbolic link: {link_path}")

    def pathtofiles(self):
        return self.temp_folder

    def print_folders(self):
        for file in natsorted(os.listdir(os.path.abspath(self.temp_folder))):
            if file.endswith(".bam"):
                print(os.path.join(self.temp_folder, file))

    def create_symbolic_links(self, file_set, subset_size):
        # Create a temporary folder
        temp_folder = tempfile.mkdtemp()

        try:
            # Ensure subset_size is not greater than the total number of files
            subset_size = min(subset_size, len(file_set))

            # Select a subset of files
            subset_files = file_set[:subset_size]

            # Create symbolic links in the temporary folder
            for file_path in subset_files:
                file_name = os.path.basename(file_path)
                link_path = os.path.join(temp_folder, file_name)

                # Create symbolic link
                os.symlink(file_path, link_path)
                print(f"Created symbolic link: {link_path}")

            print(f"\nSubset of {subset_size} symbolic links created in {temp_folder}")

        except Exception as e:
            print(f"An error occurred: {e}")

        finally:
            # Cleanup: Remove temporary folder and symbolic links in case of any exception
            shutil.rmtree(temp_folder, ignore_errors=True)
