""" Program to check whether the modified files are properly documented.
"""
import argparse
import sys
from pathlib import Path
from git_evaluation_tools import get_diff_as_json

class DoxygenMessage:
    """
    A class which stores and formats the message outputted by doxygen.
    """
    def __init__(self, file, line, level, message):
        self.file = Path(file.strip())
        self.line = int(line)
        self.level = level.strip()
        self.message = message.strip()

    def __str__(self):
        return f"{self.file} : {self.line} : {self.message}"

def unpack_doxygen_messages(doxygen_output):
    """
    Unpack the Doxygen messages from the file.

    Parameters
    ----------
    doxygen_output : str
        The file where the doxygen stderr output was printed.

    Returns
    -------
    list of DoxygenMessage
        A list of all the messages raised by doxygen.
    list of str
        A list of all the warnings raised by doxygen which are not specific to a file.
    list of str
        A list of all the errors raised by doxygen which are not specific to a file.
    """
    with open(doxygen_output, encoding="utf-8") as f:
        lines = f.readlines()
    general_errors = [l.split(':',1)[1].strip() for l in lines if l.startswith('error:')]
    general_warnings = [l.split(':',1)[1].strip() for l in lines if l.startswith('warning:')]
    messages = [DoxygenMessage(*l.split(':',3)) for l in lines if l.count(':') >= 3]
    return messages, general_warnings, general_errors

def filter_messages(doxygen_output, diff):
    """
    Filter out the Doxygen messages unrelated to this pull request.

    Parameters
    ----------
    doxygen_output : list of DoxygenMessage
        A list of all the messages raised by doxygen.

    diff : dict
        Dictionary describing all changes to the files in this pull request.

    Returns
    -------
    list of DoxygenMessage
        A list of all the messages raised by doxygen which are pertinent to this pull request.
    """
    return [o for o in doxygen_output if o.file in diff]

parser = argparse.ArgumentParser(description='Check that all new code is documented')
parser.add_argument('diff_file', metavar='diff_file', type=str,
                        help='File containing the git diff output')
parser.add_argument('doxygen_output', metavar='doxygen_output', type=str,
                        help='Stderr output from doxygen')

args = parser.parse_args()

diff = get_diff_as_json(args.diff_file)
file_messages, warnings, errors = unpack_doxygen_messages(args.doxygen_output)

filtered_messages = filter_messages(file_messages, diff)
for m in filtered_messages:
    print(m)
for w in warnings:
    print(w)
for e in errors:
    print(e)

print("Found ", len(filtered_messages)+len(warnings)+len(errors), " messages")

if filtered_messages or warnings or errors:
    sys.exit(1)
