""" Script Copies specified file types [".md", ".png", ".jpg"] from `src_dir` to `dest_dir`, preserving directory structure
    and excluding given directories.
"""
from pathlib import Path
import shutil
import argparse
from readme_to_mkdocs import replace_math_tags_with_mkdoc_compatible_tags


def copy_files(src_dir, dest_dir, exclude_dirs, extensions=(".md", ".png", ".jpg")):
    """
    Copies specified file types [".md", ".png", ".jpg"] from `src_dir` to `dest_dir`, preserving directory structure
    and excluding given directories.

    Parameters:
    - src_dir (Path): Source directory.
    - dest_dir (Path): Destination directory.
    - exclude_dirs (list): Directories to skip.
    - extensions (tuple, optional): File types to copy (default: .md, .png, .jpg).

    - Markdown files are processed with `replace_math_tags_with_mkdoc_compatible_tags`.
    - Other files are copied directly.

    Example:
    ```sh
    python docs/prepare_mkdocs.py . docs/mkdoc docs/mkdoc
    ```
    """

    exclude_dirs = [Path(e) for e in exclude_dirs]

    for e in extensions:
        for src_file in src_dir.rglob(f"*{e}"):
            if any(e in src_file.parents for e in exclude_dirs):
                print(f"Skipping: {src_file} (excluded)")
                continue

            dest_file = dest_dir / src_file.relative_to(src_dir)
            dest_file.parent.mkdir(parents=True, exist_ok=True)

            if src_file.suffix == '.md':
                replace_math_tags_with_mkdoc_compatible_tags(str(src_file), str(dest_file))
            else:
                shutil.copy2(str(src_file), str(dest_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy Markdown and PNG files while preserving directory structure.")
    parser.add_argument("src_dir", type=Path, help="Source directory")
    parser.add_argument("dest_dir", type=Path, help="Destination directory")
    parser.add_argument("exclude_dirs", nargs='*', help="Directories to exclude")
    args = parser.parse_args()

    copy_files(args.src_dir, args.dest_dir, args.exclude_dirs)
