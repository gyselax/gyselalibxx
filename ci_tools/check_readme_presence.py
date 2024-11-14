""" Program to check whether a README.md is present in any folder where new files have been added.
"""
import argparse
import glob
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from git_evaluation_tools import get_diff_as_json

gyselalibxx_root = Path(__file__).parent.parent

parser = argparse.ArgumentParser(description='Check that a README.md is present in any folder where new files have been added')
parser.add_argument('diff_file', metavar='diff_file', type=str,
                        help='File containing the git diff output')
args = parser.parse_args()

diff = get_diff_as_json(args.diff_file)

missing_readme = False

modified_folders = set(diff_file.parents[0] for diff_file in diff)

for folder in modified_folders:
    if 'src' not in folder.parts:
        continue
    expected_readme = folder.joinpath('README.md')
    if not os.path.exists(expected_readme):
        missing_readme = True
        print("Cannot find README.md file : ", expected_readme)

grep = shutil.which('grep')

for code_folder in ('src', 'simulations', 'tests', 'vendor/sll'):
    all_readme_files = glob.glob(f'{code_folder}/**/README.md', recursive=True)

    p = subprocess.run([grep, 'README.md', '-r', str(gyselalibxx_root / code_folder)], capture_output=True, check=True, encoding='utf-8')
    references = {gyselalibxx_root / readme: False for readme in all_readme_files}
    references[gyselalibxx_root / code_folder / 'README.md'] = True

    url_search = re.compile(r'\[[a-zA-Z0-9\\\_\- ]*\]\([a-zA-Z0-9\./\_\-]*README.md\)')
    for r in p.stdout.split('\n'):
        if r == '':
            continue
        calling_readme, reference_line = r.split(':',1)
        calling_folder = Path(calling_readme).parent

        for url_match in url_search.findall(reference_line):
            url = url_match.rsplit('(',1)[1].split(')')[0]
            file = (calling_folder / url).resolve()
            references[file] = True

    for readme, success in references.items():
        if not success:
            missing_readme = True
            print(f"The README.md file {readme} does not seem to be referenced from a parent folder. Please add the link to keep the website tidy.")

if missing_readme:
    sys.exit(1)
