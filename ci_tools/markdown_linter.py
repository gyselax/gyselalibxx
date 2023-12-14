"""
Script for linting .md files. This is especially important for files containing equations to ensure that the syntax can be correctly parsed by GitLab.
"""
from argparse import ArgumentParser
from collections import namedtuple
import re
import sys

import numpy as np

LineInfo = namedtuple('LineInfo', ('line', 'char'))
Expression = namedtuple('Expression', ('start_pos', 'end_pos', 'contents','closing_tag'), defaults=(None,))

def index_to_line_info(idx, text_block):
    """
    Get the line information for an index in a text block.
    """
    line = text_block[:idx].count('\n')+1
    char = idx - max(text_block[:idx].rfind('\n'), 0)+1
    return LineInfo(line, char)

# '$'
inline_math_tag = re.compile(r'(?<!\`)(?<!\$)\$(?!\`)(?!\$)')
start_inline_math_tag = re.compile(r'(?<!\$)\$\`(?!\$)')
end_inline_math_tag = re.compile(r'(?<!\$)\`\$(?!\$)')
# '$$'
standalone_math_tag = re.compile(r'(?<!\$)\$\$(?!\$)')
# '`'
inline_code_tag = re.compile(r'(?<!\`)\`(?!\`)')
multiline_code_tag = re.compile(r'\`\`\`')
# [..](..)
start_link_tag = re.compile(r'(\[)[^\[\]]*\]\([^\)]*\)')
end_link_tag = re.compile(r'\]\([^\)]*\)')

markdown_special_chars = ('*', '_', '<', '>')

unquoted_underscore = re.compile(r'(?<!\\)_')
problematic_math_underscore = re.compile(r'(?<=[\)\]}])_')
problematic_markdown_special_chars = re.compile(r'[\<\>]')

if __name__ == '__main__':
    parser = ArgumentParser(description='Check that markdown file is compatible with GitLab restrictions.')

    parser.add_argument('input_files', type=str, nargs='+', help='The file to check.')
    parser.add_argument('--fix', action='store_true', help='If this flag is used then the script will attempt to fix the file in place. Please check the output before committing')
    args = parser.parse_args()

    success = True
    warnings = False

    for file in args.input_files:
        print(f"Checking {file}")
        with open(file, 'r', encoding='utf-8') as f:
            file_contents = f.read()

        n_chars = len(file_contents)

        start_tags = [inline_math_tag, start_inline_math_tag, standalone_math_tag, inline_code_tag, multiline_code_tag, start_link_tag]
        end_tags = [inline_math_tag, end_inline_math_tag, standalone_math_tag, inline_code_tag, multiline_code_tag, end_link_tag]
        n_tag_types = len(start_tags)

        tag_names = ["inline math tag $", "inline math tag $`", "standalone equation tag $$", "inline code tag `", "multiline code tag ```", "link tag []()"]
        tags = [list(s_tag.finditer(file_contents)) for s_tag in start_tags]

        expressions = [[] for _ in range(n_tag_types+1)]
        suggested_expressions = []

        end_line_info = LineInfo(1,1)
        content_idx = 0

        while any(len(t)>0 for t in tags):
            key_idx = np.argmin([n_chars if len(t)==0 else t[0].start() for t in tags])
            start_tag = tags[key_idx][0]
            if start_tag.groups():
                found_start = start_tag.end(1)
            else:
                found_start = start_tag.end()
            start_line_info = index_to_line_info(start_tag.start(), file_contents)
            expressions[-1].append(Expression(end_line_info, start_line_info, file_contents[content_idx:start_tag.start()]))
            if key_idx in (1, 5):
                end_tag = end_tags[key_idx].search(file_contents, found_start)
            else:
                try:
                    end_tag = tags[key_idx][1]
                except IndexError:
                    end_tag = None
            if end_tag is None:
                success = False
                line_info = start_line_info
                print(f"The {tag_names[key_idx]} at line {line_info.line}, position {line_info.char} does not have a matching closing tag.", file=sys.stderr)
                tags = [[]]
            else:
                expr = file_contents[found_start:end_tag.start()]
                end_line_info = index_to_line_info(end_tag.start(), file_contents)
                expressions[key_idx].append(Expression(start_line_info, end_line_info, expr, end_tag[0]))
                content_idx = end_tag.end()
                tags = [list(s_tag.finditer(file_contents, content_idx)) for s_tag in start_tags]

        expressions[-1].append(Expression(end_line_info, index_to_line_info(n_chars, file_contents), file_contents[content_idx:]))

        content_lines = file_contents.split('\n')
        n_math_blocks = len(expressions[0]) + len(expressions[1]) + len(expressions[2])

        # Examine inline maths:
        for expr in expressions[0]:
            suggestion = expr.contents
            start_code = '$'
            end_code = '$'
            if expr.contents[0] == ' ' or expr.contents[-1] == ' ':
                success = False
                print(f"Simple inline maths expressions should not start or end with a space. ({file}: Line {expr.start_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                suggestion = suggestion.strip()
            if any(c in markdown_special_chars for c in expr.contents):
                success = False
                start_code = '$`'
                end_code = '`$'
                print(f"Simple inline maths expressions should not contain Markdown special characters {markdown_special_chars}. Please avoid these characters or use $` tags. ({file}: Line {expr.start_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
            if expr.start_pos.line != expr.end_pos.line:
                success = False
                print(f"Inline maths expressions should be written in one line. ({file}: Line {expr.start_pos.line}-{expr.end_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                suggestion = suggestion.replace('\n','')

            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, start_code+suggestion+end_code))

        # Examine inline quoted maths
        for expr in expressions[1]:
            suggestion = expr.contents
            start_code = '$`'
            end_code = '`$'
            if expr.start_pos.line != expr.end_pos.line:
                success = False
                print(f"Inline maths expressions should be written in one line. ({file}: Line {expr.start_pos.line}-{expr.end_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                suggestion = suggestion.replace('\n','')

            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, start_code+suggestion+end_code))

        # Examine standalone equations
        for expr in expressions[2]:
            suggestion = expr.contents
            start_code = '$$'
            end_code = '$$'
            if expr.start_pos.line != expr.end_pos.line and (suggestion[0]!='\n' or suggestion[-1]!='\n'):
                success = False
                print(f"Muli-line maths expressions written with one line syntax. ({file}: Line {expr.start_pos.line}-{expr.end_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                if suggestion[0]!='\n':
                    suggestion = '\n'+suggestion
                if suggestion[-1]!='\n':
                    suggestion = suggestion+'\n'

            if expr.start_pos.char != 2:
                success = False
                print(f"One-line maths expressions should be written in their own line and should start at the start of that line. ({file}: Line {expr.start_pos.line}, position {expr.start_pos.char})", file=sys.stderr)

            if expr.end_pos.char != len(content_lines[expr.end_pos.line-1]):
                success = False
                print(f"Trailing characters after a one-line maths expression. ({file}: Line {expr.start_pos.line}, position {expr.start_pos.char})", file=sys.stderr)

            if (problematic_math_underscore.search(suggestion) or problematic_markdown_special_chars.search(suggestion)) and \
                    (suggestion[0]!='\n' or suggestion[-1]!='\n'):
                success = False
                print(f"GitLab error found. Please use a multi-line maths block. ({file}: Line {expr.start_pos.line}-{expr.end_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                if suggestion[0]!='\n':
                    suggestion = '\n'+suggestion
                if suggestion[-1]!='\n':
                    suggestion = suggestion+'\n'

            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, start_code+suggestion+end_code))

        # Examine inline code tags
        for expr in expressions[3]:
            suggestion = expr.contents
            start_code = '`'
            end_code = '`'
            if expr.start_pos.line != expr.end_pos.line:
                success = False
                print(f"Inline code should be written in one line. ({file}: Line {expr.start_pos.line}-{expr.end_pos.line}, position {expr.start_pos.char})", file=sys.stderr)
                suggestion = suggestion.replace('\n','')

            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, start_code+suggestion+end_code))

        # Examine code blocks (nothing to do)
        suggested_expressions.extend(Expression(expr.start_pos, expr.end_pos, f'```{expr.contents}```') for expr in expressions[4])

        # Examine link text
        for expr in expressions[5]:
            suggestion = expr.contents
            unquoted_underscore_problems = list(unquoted_underscore.finditer(suggestion))
            if unquoted_underscore_problems:
                success = False
                for p in unquoted_underscore_problems:
                    internal_line_info = index_to_line_info(p.start(), suggestion)
                    if internal_line_info.line == 1:
                        line_info = LineInfo(expr.start_pos.line, expr.start_pos.char + internal_line_info.char-1)
                    else:
                        line_info = LineInfo(expr.start_pos.line + internal_line_info.line-1, internal_line_info.char)
                    print(f"Found non-escaped underscore. Please escape underscores and use *a* to indicate italics. ({file}: Line {line_info.line}, position {line_info.char})", file=sys.stderr)
                suggestion = unquoted_underscore.sub(r'\_', suggestion)
            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, '['+suggestion+expr.closing_tag))

        # Examine text
        for expr in expressions[6]:
            suggestion = expr.contents
            unquoted_underscore_problems = list(unquoted_underscore.finditer(suggestion))
            if unquoted_underscore_problems:
                success = False
                for p in unquoted_underscore_problems:
                    internal_line_info = index_to_line_info(p.start(), suggestion)
                    if internal_line_info.line == 1:
                        line_info = LineInfo(expr.start_pos.line, expr.start_pos.char + internal_line_info.char-1)
                    else:
                        line_info = LineInfo(expr.start_pos.line + internal_line_info.line-1, internal_line_info.char)
                    print(f"Found non-escaped underscore. Please escape underscores and use *a* to indicate italics. ({file}: Line {line_info.line}, position {line_info.char})", file=sys.stderr)
                suggestion = unquoted_underscore.sub(r'\_', suggestion)
            suggested_expressions.append(Expression(expr.start_pos, expr.end_pos, suggestion))

        suggested_expressions.sort(key=lambda expr: (expr.start_pos.line, expr.start_pos.char))

        if n_math_blocks > 50:
            print(f"File contains {n_math_blocks} math blocks but GitLab will only render 50. Some equations may not show up correctly on GitLab", file=sys.stderr)
            warnings = True

        if args.fix:
            print(f"Fixing {file}")
            out_code = ''.join(e.contents for e in suggested_expressions)
            with open(file, 'w', encoding='utf-8') as f:
                print(out_code, end='', file=f)

    if not success:
        sys.exit(1)
    elif warnings:
        sys.exit(2)
    else:
        sys.exit(0)
