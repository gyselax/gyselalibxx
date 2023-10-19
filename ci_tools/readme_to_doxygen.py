""" Script for translating Markdown files compatible with GitHub/GitLab to a format that can be understood by Doxygen.
"""
from argparse import ArgumentParser
import os
import re

# '$', not preceded by '$' or followed by '$'
single_math_tag = re.compile(r'(?<!\$)\$(?!\$)')
# '$$', not preceded by '$' or followed by '$'
multiline_math_tag = re.compile(r'(?<!\$)\$\$(?!\$)')
reference_tag = re.compile(r'(\[[^\]]*\])(\([^\)]*\))')
section_tag = re.compile(r'\n## ([^\n]*)\n')

doxygen_tag_dict = {" ":"_", "+":"x", os.path.sep: '_'}

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

def format_equations(line, start_tag, end_tag, start_replace, end_replace, contents = None, idx = None):
    """
    Format equations with Doxygen style.

    Use regexes to find equations. If they are not in code blocks (e.g. for
    documentation explaining how to format equations) then replace Markdown
    formatting with Doxygen style formatting.

    Parameters
    ----------
    line : str
        A line of markdown text.

    start_tag : re.Pattern
        A regex pattern matching the start of the equation.

    end_tag : re.Pattern
        A regex pattern matching the end of the equation.

    start_replace : str
        A Doxygen-compatible string which should be used to introduce the equation.

    end_replace : str
        A Doxygen-compatible string which should be used to close the equation.

    contents : list, optional
        A list of all the lines in the file. This should be provided for multi-line
        equations.

    idx : int
        The index of the line in the contents. This should be provided for multi-line
        equations.

    Returns
    -------
    str
        The updated line.

    int
        The updated line index.
    """
    code_blocks = get_code_blocks(line)
    start_match = start_tag.search(line)
    if contents:
        n = len(contents)

    while start_match:
        if not any(start < start_match.start() < end for start, end in code_blocks):
            end_match = end_tag.search(line, start_match.end())
            if end_match is None:
                assert contents is not None

            while end_match is None and idx < n-1:
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

root_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))

if __name__ == '__main__':
    parser = ArgumentParser(description='Tool for translating Markdown files to a format that can be understood by Doxygen.')

    parser.add_argument('input_file', type=str, help='The file to treated.')
    parser.add_argument('output_file', type=str, help='The file where the results should be saved.')
    args = parser.parse_args()

    file = args.input_file
    folder = os.path.dirname(file)

    with open(file, 'r', encoding='utf-8') as f:
        contents = f.readlines()

    title = contents[0]

    assert title.startswith('#')

    title = title[1:].strip()

    relpath = os.path.relpath(file, start=root_dir)

    if file.endswith('README.md'):
        reldir = os.path.dirname(relpath)
        tag = get_compatible_tag(reldir)
    else:
        tag = get_compatible_tag(relpath[:-3])

    # Add tag to page
    if tag:
        title = f'@page {tag} {title}\n'
    else:
        title = '@mainpage\n'

    body = []
    in_code = False
    in_math_code = False
    line_index = 1
    n_lines = len(contents)
    while line_index < n_lines:
        line = contents[line_index]
        stripped = line.strip()
        if stripped.startswith('```'):
            if in_code:
                in_code = False
            elif in_math_code:
                body.append('@f]')
                in_math_code = False
                continue
            elif stripped == '```math':
                body.append('@f[')
                in_math_code = True
                continue
            else:
                in_code = True
        elif not in_code:
            if line.startswith('##') and line[2] != '#':
                sec_title = line[2:].strip()
                sec_tag = get_compatible_tag(sec_title)
                line = f"\n@section {sec_tag} {sec_title}\n"
            else:
                # Replace inline math with Doxygen tag
                line, _ = format_equations(line, single_math_tag, single_math_tag, "@f$", "@f$")
                line, line_index = format_equations(line, multiline_math_tag, multiline_math_tag, "@f[", "@f]", contents, line_index)

                # Look for references
                match_found = reference_tag.search(line)
                while match_found:
                    # Get path to referenced file
                    path = match_found.group(2)[1:-1]

                    # Replace references to READMEs with the page tag
                    if path.endswith('README.md') and not os.path.isabs(path):
                        # Get the path to the referenced file directory from the root directory
                        path = os.path.normpath(os.path.join(folder, path))
                        relpath = os.path.relpath(path, start=root_dir)
                        reldir = os.path.dirname(relpath)

                        # Get the path to the referenced file directory from the file containing the reference
                        path_from_file = os.path.relpath(os.path.dirname(path), start=os.path.dirname(file))

                        # Check if the referenced file is a subpage.
                        # It is a subpage if it is found in a sub-folder and no other README is found in a
                        # closer sub-folder
                        is_subpage = '..' not in path_from_file
                        if is_subpage:
                            intermediate_folders = path_from_file.split('/')
                            for i in range(1,len(intermediate_folders)-1):
                                # Check if a README can be found in the intermediate file
                                if os.path.exists(os.path.join(folder, *intermediate_folders[:i], "README.md")):
                                    is_subpage = False
                                    break

                        # Calculate the tag of the page
                        subpage_tag = get_compatible_tag(reldir)

                        # Create the complete doxygen command to reference the page.
                        doxygen_command = "@subpage " if is_subpage else "@ref "
                        doxygen_ref = doxygen_command + subpage_tag

                        # Replace the markdown command with the doxygen command
                        line = line[:match_found.start()] + doxygen_ref + line[match_found.end():]
                        match_found = reference_tag.search(line, match_found.start() + len(doxygen_ref))
                    else:
                        match_found = reference_tag.search(line, match_found.end())

        body.append(line)
        line_index += 1


    directory = os.path.dirname(args.output_file)
    os.makedirs(directory, exist_ok = True)

    with open(args.output_file, 'w', encoding='utf-8') as f:
        f.write(title)
        f.write(''.join(body))
