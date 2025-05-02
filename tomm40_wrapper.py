#!/usr/bin/env python
import os
import sys
import subprocess

def main():
    # Get the environment prefix (root of the venv)
    env_root = sys.prefix

    # Construct the path to the actual script and data directory
    # data_files in setup.py installs relative to sys.prefix
    data_dir = os.path.join(env_root, "opt", "tomm40_wgs")
    real_script_path = os.path.join(data_dir, "TOMM40_WGS")

    if not os.path.exists(real_script_path):
        # Let's check an alternative common location just in case
        alt_data_dir = os.path.join(env_root, "lib", f"python{sys.version_info.major}.{sys.version_info.minor}", "opt", "tomm40_wgs")
        alt_real_script_path = os.path.join(alt_data_dir, "TOMM40_WGS")
        if os.path.exists(alt_real_script_path):
            real_script_path = alt_real_script_path
            data_dir = alt_data_dir
        else:
            print(f"Error: Could not find the main TOMM40_WGS script.", file=sys.stderr)
            print(f"       Checked: {real_script_path}", file=sys.stderr)
            print(f"       Checked: {alt_real_script_path}", file=sys.stderr)
            sys.exit(1)

    # Set an environment variable for the real script to use
    env = os.environ.copy()
    env["TOMM40_DATA_DIR"] = data_dir

    # Get command line arguments passed to the wrapper
    args = [real_script_path] + sys.argv[1:]

    # Execute the real script
    try:
        # Use execvpe to replace the wrapper process with the real script
        # This ensures signals (like Ctrl+C) are handled correctly by the real script
        os.execvpe(real_script_path, args, env)
    except OSError as e:
        print(f"Error executing {real_script_path}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()


