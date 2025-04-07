import os
import subprocess

# Get the current working directory
base_directory = os.getcwd()

# Iterate through each folder in the base directory
for folder_name in os.listdir(base_directory):
    folder_path = os.path.join(base_directory, folder_name)
    base_folder_name = os.path.basename(folder_path)

    # Check if the item is a directory and matches the pattern
    if os.path.isdir(folder_path) and folder_name.startswith('D_LANDAU4D_FFT_'):
        # List all files in the folder
        files = os.listdir(folder_path)

        # Find the first file that matches the pattern g*.dat
        for file_name in files:
            if ( file_name.startswith('g') or file_name.startswith('c') ) and file_name.endswith('.dat'):
                # Construct the input and output file paths
                input_file_path = os.path.join(folder_path, file_name)
                output_file_name = f'timer.{folder_name}'
                output_file_path = os.path.join(base_directory, output_file_name)

                # Run the kp_reader command
                command = f'../../kokkos-tools/install/bin/kp_reader {input_file_path} > {output_file_path}'
                # Run the kp_reader command
                try:
                    subprocess.run(command, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Command failed with error: {e}")

                # Break after processing the first matching file
                break

