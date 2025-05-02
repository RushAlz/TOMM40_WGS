import setuptools
import os

# --- Configuration ---
package_name = "tomm40_wgs"
package_version = "0.0.1"
source_root_dir = "." # Assumes setup.py is in the TOMM40_WGS directory
wrapper_module_name = "tomm40_wrapper"
executable_name = "TOMM40_WGS"
original_script_name = "TOMM40_WGS" # The actual script file

# --- Helper function to find data files ---
def find_data_files(source_dir, install_base_dir):
    data_files_spec = []
    # Walk through the source directory
    for dirpath, dirnames, filenames in os.walk(source_dir):
        # Exclude hidden directories like .git if they exist
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]

        # Determine the target directory relative to install_base_dir
        relative_subdir = os.path.relpath(dirpath, source_dir)
        target_dir = os.path.join(install_base_dir, relative_subdir) if relative_subdir != "." else install_base_dir

        # List of source files in this directory
        sources = []
        for filename in filenames:
            # Exclude setup.py itself
            if dirpath == source_dir and filename == "setup.py":
                 continue
            # Exclude the wrapper script module
            if dirpath == source_dir and filename == f"{wrapper_module_name}.py":
                 continue
            # Exclude hidden files
            if filename.startswith("."):
                continue
            # Add other exclusions if needed (e.g., __pycache__)
            if "__pycache__" in dirpath:
                continue

            sources.append(os.path.join(dirpath, filename))

        if sources:
            data_files_spec.append((target_dir, sources))
    return data_files_spec

# --- Determine data files installation --- 
# Installs all files (including the original script) into opt/tomm40_wgs
data_install_dir_relative = os.path.join("opt", package_name)
package_data_files = find_data_files(source_root_dir, data_install_dir_relative)

# --- Read README for long description ---
try:
    with open(os.path.join(source_root_dir, "README.md"), "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "TOMM40 WGS analysis pipeline. See homepage for details."

# --- Setup configuration ---
setuptools.setup(
    name=package_name,
    version=package_version,
    author="Ricardo A. Vialle",
    author_email="ricardo_a_vialle@rush.edu",
    description="Genotyping TOMM40 '523 Poly-T Polymorphisms from Whole-Genome Sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RushAlz/TOMM40_WGS/",
    license="GPL-2.0-or-later",

    # Explicitly state no Python packages are being installed (except the wrapper)
    packages=[],
    # Specify the wrapper module
    py_modules=[wrapper_module_name],

    # Use entry points to create the executable command from the wrapper
    entry_points={
        'console_scripts': [
            f'{executable_name} = {wrapper_module_name}:main',
        ],
    },

    # Specify Python dependencies
    install_requires=[
        "pandas>=2.2.3",
        "snakemake>=9.3.0"
    ],

    # Specify data files to be installed (includes the original script)
    data_files=package_data_files,

    # Classifiers help users find your project
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],

    # Minimum Python version requirement
    python_requires=">=3.9",

    # Project URLs
    project_urls={
        "Bug Tracker": "https://github.com/RushAlz/TOMM40_WGS/issues",
        "Source Code": "https://github.com/RushAlz/TOMM40_WGS/",
    },

    # Include package data (alternative way, might not be needed with data_files)
    # include_package_data=True, 
    # package_data={ '': ['*'] }, # Be careful with this, data_files is more explicit

    # Note: This setup uses a wrapper script.
    # The *original* TOMM40_WGS script needs a one-time modification
    # to read the TOMM40_DATA_DIR environment variable set by the wrapper.
)

print("\n" + "*"*70)
print("IMPORTANT NOTES FOR WRAPPER SOLUTION:")
print("- This setup uses a Python wrapper script.")
print("- The original `TOMM40_WGS` script (the one in your source directory)")
print("  needs to be modified *once* before you run `pip install`.")
print("  It should be changed to read the data directory path from the")
print("  environment variable 'TOMM40_DATA_DIR' instead of using hardcoded")
print("  paths or paths relative to the script's original location.")
print("  For example, in bash, replace lines setting SCRIPT_DIR or similar with:")
print("    DATA_DIR=\"${TOMM40_DATA_DIR:-.}\" # Use env var, default to current dir if unset")
print("    SNAKEFILE=\"$DATA_DIR/Snakefile\"")
print("*"*70 + "\n")


