from pathlib import Path
import shutil
import argparse
import sys
from readme_to_mkdocs import replace_math_tags_with_mkdoc_compatible_tags

def copy_files(src_dir, dest_dir, exclude_dirs, extensions=(".md", ".png")):
    exclude_dirs = [Path(e) for e in exclude_dirs]
    for e in extensions:
        for src_file in src_dir.rglob(f"*{e}"):
            if any(e in src_file.parents for e in exclude_dirs):
                continue
            # Determine relative path and corresponding destination
            dest_file = dest_dir / src_file

            dest_file.parent.mkdir(parents=True, exist_ok=True)

            if src_file.suffix == '.md':
                replace_math_tags_with_mkdoc_compatible_tags(str(src_file), str(dest_file))
            else:
                shutil.copy2(str(src_file), str(dest_file))
            print(f"Copied: {src_file} -> {dest_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy Markdown and PNG files while preserving directory structure.")
    parser.add_argument("src_dir", type=Path, help="Source directory to search")
    parser.add_argument("dest_dir", type=Path, help="Destination directory to copy files to")
    parser.add_argument("exclude_dirs", nargs='*', help="Directories that should be ignored during the copy")
    args = parser.parse_args()
    
    copy_files(args.src_dir, args.dest_dir, args.exclude_dirs)
