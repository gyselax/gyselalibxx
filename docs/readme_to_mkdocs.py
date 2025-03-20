""" Script for translating Markdown files compatible with GitHub/GitLab to a format that can be understood by Doxygen.
"""
from argparse import ArgumentParser
import os
import re

# '$', not preceded by '$' or followed by '$'
single_math_start_tag = re.compile(r'(?<!\$)\$\`(?!\$)')
single_math_end_tag = re.compile(r'(?<!\$)\`\$(?!\$)')
single_math_tag = re.compile(r'(?<!@f)(?<!\$)\$(?!\$)')
# '$$', not preceded by '$' or followed by '$'
multiline_math_tag = re.compile(r'(?<!\$)\$\$(?!\$)')
reference_tag = re.compile(r'(\[[^\]]*\])(\([^\)]*\))')
section_tag = re.compile(r'\n## ([^\n]*)\n')

doxygen_tag_dict = {" ":"_", "+":"x", "'":'_', os.path.sep: '_', '?':'', '!':'', ',':'', r'\_': '_'}

def get_compatible_tag(tag):
    """
    Get a Doxygen tag which doesn't contain problematic symbols.

    Construct a Doxygen tag from a string by replacing any problematic
    symbols.

    Parameters
    ----------
    tag : str
        The string from which the tag should be constructed.

    Returns
    -------
    str
        A string that can be used as a Doxygen tag.
    """
    for k, v in doxygen_tag_dict.items():
        tag = tag.replace(k, v)
    return tag

def get_code_blocks(line):
    """
    Get the position of code blocks in the line.

    Parameters
    ----------
    line : str
        A line of markdown text.

    Returns
    -------
    list of 2-tuples
        Returns a list of 2-tuples describing the start and end of the code blocks in the line.
    """
    # Find the code blocks in the line
    code_tags = [i for i,l in enumerate(line) if l=='`']
    n_code_tags = len(code_tags)
    # Ensure that all inline code blocks are closed
    assert n_code_tags % 2 == 0
    n_code_tags //= 2
    return [(code_tags[2*i], code_tags[2*i+1]) for i in range(n_code_tags)]

def format_equations(line, tags, replacements, contents=None, idx=None):
    """
    Format equations with Doxygen style.

    Parameters
    ----------
    line : str
        A line of markdown text.

    tags : tuple (re.Pattern, re.Pattern)
        A tuple containing regex patterns matching the start and end of the equation.

    replacements : tuple (str, str)
        A tuple containing Doxygen-compatible strings to introduce and close the equation.

    contents : list, optional
        A list of all the lines in the file. Required for multi-line equations.

    idx : int, optional
        The index of the line in the contents. Required for multi-line equations.

    Example
    --------

    start_tag = re.compile(r"$$")
    end_tag = re.compile(r"$$")
    formatted_line, new_idx = format_equations(line, (start_tag, end_tag), ("\\(", "\\)"), contents, idx)

    Returns
    -------
    str
        The updated line.

    int
        The updated line index.
    """
    start_tag, end_tag = tags
    start_replace, end_replace = replacements

    code_blocks = get_code_blocks(line)
    start_match = start_tag.search(line)
    n = len(contents) if contents else 1

    while start_match:
        if not any(start < start_match.start() < end for start, end in code_blocks):
            end_match = end_tag.search(line, start_match.end())
            if end_match is None:
                assert contents is not None

            while end_match is None and idx < n - 1:
                idx += 1
                line += '\n' + contents[idx]
                end_match = end_tag.search(line, start_match.end())

            replacement = start_replace + line[start_match.end():end_match.start()] + end_replace
            line = line[:start_match.start()] + replacement + line[end_match.end():]
            code_blocks = get_code_blocks(line)
            start_match = start_tag.search(line, start_match.start() + len(replacement))
        else:
            start_match = start_tag.search(line, start_match.end())

    return line, idx


def replace_math_tags_with_mkdoc_compatible_tags(in_file, out_file):
    """
    Converts math block tags in a Markdown file to MkDocs-compatible tags.

    This function reads a Markdown file, detects math block tags formatted as ` ```math `,
    and replaces them with MkDocs-compatible tags. Inline math expressions
    are also converted accordingly.

    Args:
        in_file (str): Path to the input Markdown file.
        out_file (str): Path to the output file where the modified content will be written.
    """
    print(in_file, out_file)
    with open(in_file, 'r', encoding='utf-8') as f:
        contents = f.readlines()

    body = []
    in_code = False
    in_math_code = False
    for line in contents:
        stripped = line.strip()
        if stripped.startswith('```'):
            if in_code:
                in_code = False
            elif in_math_code:
                # Replace block tag math with mkdoc tag
                body.append('\\]\n')
                in_math_code = False
                continue
            elif stripped == '```math':
                # Replace block tag math with mkdoc tag
                body.append('\n\\[\n')
                in_math_code = True
                continue
            else:
                in_code = True
        elif not in_code and not in_math_code:
            # Replace inline math with mkdoc tag
            line, _ = format_equations(line, (single_math_start_tag, single_math_end_tag), ("\\(", "\\)"))

        body.append(line)

    directory = os.path.dirname(out_file)
    os.makedirs(directory, exist_ok = True)

    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(''.join(body))

if __name__ == '__main__':
    parser = ArgumentParser(description='Tool for translating Markdown files to a format that can be understood by Doxygen.')

    parser.add_argument('input_file', type=str, help='The file to treated.')
    parser.add_argument('output_file', type=str, help='The file where the results should be saved.')
    args = parser.parse_args()

    replace_math_tags_with_mkdoc_compatible_tags(args.input_file, args.output_file)
