""" Program to check whether a README.md is present in any folder where new files have been added.
"""
import argparse
import os
import sys
from git_evaluation_tools import get_diff_as_json


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

if missing_readme:
    sys.exit(1)
