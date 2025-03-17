import os
import shutil
import argparse

def copy_files(src_dir, dest_dir, extensions=(".md", ".png")):
    # Get a list of all subdirectories
    subdirs = [os.path.join(root) for root, _, _ in os.walk(src_dir)]
    subdirs.append(src_dir)
    for root in subdirs:
        
        # List all files in the current directory
        files = [f for f in os.listdir(root) if os.path.isfile(os.path.join(root, f))]
        
        # Filter files with the desired extensions
        relevant_files = [f for f in files if f.lower().endswith(extensions)]
        
        if relevant_files:
            # Determine relative path and corresponding destination
            rel_path = os.path.relpath(root, src_dir)
            target_dir = os.path.join(dest_dir, rel_path)
            os.makedirs(target_dir, exist_ok=True)
            
            for file in relevant_files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(target_dir, file)
                shutil.copy2(src_file, dest_file)
                print(f"Copied: {src_file} -> {dest_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy Markdown and PNG files while preserving directory structure.")
    parser.add_argument("src_dir", help="Source directory to search")
    parser.add_argument("dest_dir", help="Destination directory to copy files to")
    args = parser.parse_args()
    
    copy_files(args.src_dir, args.dest_dir)
