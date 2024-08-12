""" A program which uses cppcheck to ensure that gyselalibxx follows good coding practices.
The dump option is used to write specific rules for gyselalibxx.
The rules that are currently implemented are:
    - Check for forbidden DDC keywords
    - Check for unnecessary auto type
    - SCheck for function returning a FieldMem saved into a Field
    - Check for use of Field<double... or FieldMem<double... instead of DField/DFieldMem
    - Check for continuous dimensions with unnecessarily long names
    - Check for grids which are not named as such
    - Check for missing #pragma once in a .hpp file
    - Check that "" is used to include files from the Gysela library and <> is used to include files from other libraries
"""
import argparse
from itertools import chain
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np
FATAL = "\033[1m\033[31mfatal\033[0m"
ERROR = "\033[1m\033[31merror\033[0m"
STYLE = "\033[1m\033[33mstyle error\033[0m"

possible_error_levels = {STYLE: 1, ERROR: 2, FATAL: 3}

error_level = 0

ddc_keyword_map = {'DiscreteElement': 'Idx',
                   'DiscreteVector' : 'IdxStep',
                   'DiscreteDomain' : 'IdxRange',
                   'Chunk'          : 'FieldMem',
                   'ChunkSpan'      : 'Field',
                   'ChunkView'      : 'ConstField',
                   'span_view'      : 'get_field',
                   'span_cview'     : 'get_const_field',
                   'domain'         : 'get_idx_range',
                   'spline_domain'  : 'get_spline_idx_range',
                   'UniformPointSampling' : 'UniformGridBase',
                   'NonUniformPointSampling' : 'NonUniformGridBase'}

mirror_functions = {'create_mirror', 'create_mirror_and_copy', 'create_mirror_view', 'create_mirror_view_and_copy'}

HOME_DIR = Path(__file__).parent.parent.absolute()
global_folders = [HOME_DIR / f for f in ('src', 'simulations', 'tests')]

auto_functions = set()
field_mem_functions = set()

def report_error(level, file, linenr, message):
    """
    A tool to report any errors.

    Parameters
    ----------
    level : str
        One of STYLE, ERROR, FATAL.
    file : FileObject
        The file where the error is found.
    linenr : int
        The line number where the error was raised from
    message : str
        The error message.
    """
    global error_level #pylint: disable=global-statement
    error_level = max(error_level, possible_error_levels[level])
    print(f'{level} {file.file} : {linenr} : {message}')

class FileObject:
    """
    A class containing the cppcheck output which describes a file.

    A class which parses the cppcheck dump file describing a file. It unpacks useful information
    from this, such as the variable descriptors and a list splitting the code into lines.

    Parameters
    ----------
    file : pathlib.Path
        A path to the dump file.
    """
    def __init__(self, file):
        self._file = file.parent / file.stem
        tree = ET.parse(file)
        self.root = tree.getroot()
        self.raw = [child.attrib for child in self.root.findall('./rawtokens/tok')]
        self.data = [child.attrib for child in self.root.findall('./dump/tokenlist/token')]
        self.raw_xml = self.root.find('./rawtokens')
        self.data_xml = self.root.find('./dump/tokenlist')
        self.variables = [child.attrib for child in self.root.findall('./dump/variables/var')]

        code_lines = []
        current_line = []
        for i,c in enumerate(self.data):
            if c['str'] in ('{', '}', ';'):
                code_lines.append(current_line)
                current_line = []
            else:
                current_line.append(c)

        self.lines = code_lines

        using_indices = [i for i, d in enumerate(self.raw) if d['str'] == 'using']
        expr_end = [i for i, d in enumerate(self.raw) if d['str'] == ';']
        self.aliases = {self.raw[idx+1]['str']: self.raw[idx+3:expr_end[np.searchsorted(expr_end, idx)]]
                    for idx in using_indices}

        self._identify_potential_bugs()

    @property
    def file(self):
        """
        The path to the file being analysed.
        """
        return self._file

    def _identify_potential_bugs(self):
        """
        A function called by the constructor which identifies functions which may cause bugs or false positives.
        This means searching for functions which return a FieldMem or an auto. If the function returns auto
        then there may be a false positive for an unnecessary auto when it is called. If the function returns a
        FieldMem then it may cause bugs if it is not saved into the same type.
        """
        for line in self.lines:
            func_idx = next((i for i,v in enumerate(line) if 'function' in v), None)
            if func_idx is not None:
                func_name = line[func_idx]['str']
                return_type = ''.join(v['str'] for v in line[:func_idx])
                if line[func_idx-1]['str'] == 'auto':
                    auto_functions.add(func_name)
                elif re.search(r'\b(ddc::Chunk)|(FieldMem)\b', return_type):
                    field_mem_functions.add(func_name)

def get_relevant_files():
    """
    Search through the project to find and group files to be analysed. The files are grouped by geometry to
    prevent errors about doubly defined aliases between geometries.
    """
    src_geometries = {f.stem: f.absolute() for f in (HOME_DIR / 'src').glob('geometry*')}
    simu_geometries = {f.stem: f.absolute() for f in (HOME_DIR / 'simulations').glob('geometry*')}
    test_geometries = {f.stem: f.absolute() for f in (HOME_DIR / 'tests').glob('geometry*')}
    assert all(g in src_geometries for g in simu_geometries)
    assert all(g in src_geometries for g in test_geometries)
    test_geometries['multipatch'] = HOME_DIR / 'tests/multipatch/'

    possible_geometries = {'src': src_geometries,
                           'simulations': simu_geometries,
                           'tests': test_geometries}

    relevant_files = {}
    # Collect all files in tests/ and simulations/ and group them by geometry
    for top_level_folder in ('tests', 'simulations'):
        files = set((HOME_DIR / top_level_folder).glob('**/*.hpp')).union((HOME_DIR / top_level_folder).glob('**/*.cpp'))

        for g, folder in possible_geometries[top_level_folder].items():
            relevant_files_in_folder = set(f for f in files if folder in f.parents)
            files = files - relevant_files_in_folder
            relevant_files.setdefault(g, set()).update(relevant_files_in_folder)

        # Files not associated with a particular geometry are saved as 'general'
        relevant_files.setdefault('general', set()).update(files)

    files = set((HOME_DIR / 'src').glob('**/*.hpp')).union((HOME_DIR / 'src').glob('**/*.cpp'))

    # Collect all files in src/ and group them by geometry
    for g, folder in possible_geometries['src'].items():
        relevant_files_in_folder = set(f for f in files if folder in f.parents)
        files = files - relevant_files_in_folder
        relevant_files.setdefault(g, set()).update(relevant_files_in_folder)

    # Any remaining files from src not associated with a particular geometry are relevant to all configurations
    for v in relevant_files.values():
        v.update(files)

    return relevant_files

def check_for_ddc_keywords(file):
    """
    Test to ensure DDC keywords are not used.

    This test excludes certain files which cannot use ddc_aliases. These files are:
    - ddc_aliases.hpp
    - ddc_helper.hpp
    - transpose.hpp
    - src/data_types/*.hpp

    The forbidden keywords are defined in the dictionary ddc_keyword_map at the top of this file.
    """
    file_exceptions = [HOME_DIR / 'src' / 'utils' / 'ddc_helper.hpp',
                       HOME_DIR / 'src' / 'utils' / 'ddc_aliases.hpp',
                       HOME_DIR / 'src' / 'utils' / 'transpose.hpp']
    if (HOME_DIR / 'src' / 'data_types') in file.file.parents or file.file in file_exceptions:
        return

    ddc_usage = [d.attrib for k in ddc_keyword_map for d in file.raw_xml.findall(".tok[@str='"+k+"']")]
    ddc_usage = [d for d in ddc_usage if d['str'] not in ('domain', 'spline_domain') \
                                        or file.raw[file.raw.index(d)-1]['str'] == '.']

    for d in ddc_usage:
        key = d['str']
        new_key = ddc_keyword_map[key]
        report_error(STYLE, file, d['linenr'], f"Please use the new naming conventions instead of DDC names ({key}->{new_key})")

def search_for_unnecessary_auto(file):
    """
    Test to ensure variable type is not 'auto' unless saving a function returning auto.
    """
    # Find all variables in the file with a type of 'auto'
    variables = [v for v in file.variables if v['typeStartToken'] == v['typeEndToken']]
    variable_types = [file.data_xml.find(".token[@id='"+v['typeStartToken']+"']") for v in variables]
    auto_variables = [file.data_xml.find(".token[@id='"+v['nameToken']+"']") \
                      for v,t in zip(variables, variable_types) if t is not None and t.attrib['str'] == 'auto']

    # Find the name and location of the first use of the auto variables
    var_data_idx = [file.data.index(v.attrib) for v in auto_variables if v]
    var_names = [v.attrib['str'] for v in auto_variables if v]
    for idx, var_name in enumerate(var_names):
        start = var_data_idx[idx]
        end = next(i for i,v in enumerate(file.data[start:], start) if v['str'] == ';')
        # If declaration is of the form auto x = ... instead of auto x(...) then examine next code phrase
        if end == start + 1:
            assert file.data[end+1]['str'] == var_name
            start = end + 1
            end = next(i for i,v in enumerate(file.data[start:], start) if v['str'] == ';')

        # If rhs of assignment is not a mirror function (where auto is compulsory) or a function returning auto
        # then raise an error
        if not any(v['str'] in chain(mirror_functions, auto_functions) for v in file.data[start:end]):
            report_error(ERROR, file, file.data[start]['linenr'], f"Please use explicit types instead of auto ({var_name})")

def search_for_bad_create_mirror(file):
    """
    Test for instances where create_mirror functions are called incorrectly.

    An incorrect create_mirror call is one which does not save the result to an auto type (leading to CPU/GPU problems) or
    one which passes Kokkos::DefaultHostExecutionSpace() as the first argument (This can cause synchronicity errors leading
    to wrong results).
    """
    for line in file.lines:
        if any(v['str'] in mirror_functions for v in line):
            first_elem = line[0]
            if 'variable' not in first_elem and first_elem['str'] == 'auto':
                report_error(FATAL, file, first_elem['linenr'], "Mirror function should be saved into auto type to handle CPU and GPU.")
            open_bracket = next(i for i,v in enumerate(line) if v['str'] == '(')
            comma = next((i for i,v in enumerate(line[open_bracket:], open_bracket) if v['str'] == ','), len(line)-1)
            first_arg = ''.join(v['str'] for v in line[open_bracket+1:comma])
            if first_arg == 'Kokkos::DefaultHostExecutionSpace()':
                line_str = ' '.join(l['str'] for l in line)
                msg = ('Possible synchronicity error. Kokkos::DefaultHostExecutionSpace() should never be '
                      f'the first argument to a mirror function ({line_str})')
                report_error(FATAL, file, first_elem['linenr'], msg)

def search_for_bad_memory(file):
    """
    Test to find dangling memory pointers.

    This test looks for places where one of the functions returning a FieldMem is used. If the result is not saved into
    a FieldMem type then the data will be deallocated leading to a memory error somewhere much later in the program when
    the Field is used.
    """
    for line in file.lines:
        func_idx = next((i for i,v in enumerate(line) if v['str'] in field_mem_functions), None)
        if func_idx is None:
            continue
        func_use_id = line[func_idx]['id']
        if file.root.find(f"./dump/scopes//function[@tokenDef='{func_use_id}']"):
            # Function definition
            continue
        first_var_idx = next((i for i,v in enumerate(line) if 'variable' in v), None)
        if first_var_idx is None or line[first_var_idx+1]['str'] not in ('(', '='):
            # No variable found
            continue
        first_var = line[first_var_idx]
        linenr = first_var['linenr']
        scope = first_var['scope']
        var_id = first_var['variable']
        var_tag = file.root.find(f".//var[@scope='{scope}'][@id='{var_id}']")
        while var_tag is None:
            scope = file.root.find(f".//scopes/scope[@id='{scope}']").attrib['nestedIn']
            var_tag = file.root.find(f".//var[@scope='{scope}'][@id='{var_id}']")

        var_description = var_tag.attrib
        start_idx = file.data.index(file.root.find(f".//token[@id='{var_description['typeStartToken']}'][@scope='{scope}']").attrib)
        end_idx = file.data.index(file.root.find(f".//token[@id='{var_description['typeEndToken']}'][@scope='{scope}']").attrib)
        type_name = ''.join(v['str'] for v in file.data[start_idx:end_idx+1])
        if 'Mem' not in type_name:
            line_str = ' '.join(l['str'] for l in line)
            msg = ( 'Possible memory error. Mem not found in type where a FieldMem object is saved to.',
                    'Returned FieldMem object must be saved in a FieldMem else data may be unintentionally deallocated.',
                    f'({line_str})')
            report_error(FATAL, file, linenr, msg)

def search_for_bad_aliases(file):
    """
    Test to find instances of Field<double/FieldMem<double and recommend the use of DField/DFieldMem.
    """
    forbidden_substrings = ('Chunk', 'Span', 'DDom', 'IDom')
    for a_name, val in file.aliases.items():
        val_str = ''.join(v['str'] for v in val)
        if re.search(r'\bField<(const )?double', val_str) or \
                re.search(r'\bFieldMem<(const )?double', val_str):
            report_error(STYLE, file, val[0]['linenr'], f'double does not need to be a template parameter. Please use DField/DFieldMem ({val_str})')
        if any(bad_name in a_name for bad_name in forbidden_substrings):
            report_error(STYLE, file, val[0]['linenr'], f"Name seems to contain DDC keywords. Please use the new naming conventions (using {a_name} = {val_str})")

    template_declarations = file.root.findall(".//TokenAndName")
    for temp in template_declarations:
        if Path(temp.attrib['file']) != file.file:
            continue
        possible_nodes = list(temp)
        template_class_loc = [possible_nodes.index(tok) \
                              for tok in temp.findall("./template-token[@str='class']")]
        template_class_loc = [l+2 if possible_nodes[l+1].attrib['str'] == '...' else l+1 \
                              for l in template_class_loc]
        template_class_names = [possible_nodes[l].attrib['str'] for l in template_class_loc]

        for name in template_class_names:
            if any(bad_name in name for bad_name in forbidden_substrings):
                msg = f"Name seems to contain DDC keywords. Please use the new naming conventions ({name})"
                report_error(STYLE, file, temp.attrib['line'], msg)
            if name == 'Grid':
                msg = "The name 'Grid' is too general and should be avoided as it implies that the type can be multi-D. Please use Grid1D"
                report_error(STYLE, file, temp.attrib['line'], msg)

def check_struct_names(file):
    """
    Test the names of structures. If they do not follow the new coding conventions then a style error is raised.

    The error is raised if:
    - a dimension is called DimX or RDimX (ie starts with Dim or RDim followed by a non-numeric name.
    - a grid does not have Grid in the name
    """
    for line in file.lines:
        if line and line[0]['str'] == 'struct':
            struct_name = line[1]['str']
            if (struct_name.startswith('Dim') and struct_name[3:].isalpha()) or \
                    (struct_name.startswith('RDim') and struct_name[4:].isalpha()):
                report_error(STYLE, file, line[0]['linenr'], f'Dim prefix is unnecessary for struct ({struct_name})')
            if len(line) > 3 and line[2] == ':' and any('GridBase' in v['str'] for v in line) and 'Grid' not in struct_name:
                report_error(STYLE, file, line[0]['linenr'], f'Grid missing from struct name ({struct_name})')

def search_for_uid(file):
    """
    Test for use of uid(). The DDC internal unique identifier should never be used directly.
    """
    if not file.data_xml:
        return
    uid_usage = [d.attrib for d in file.data_xml.findall(".token[@str='uid']")]
    uid_indices = [file.data.index(d) for d in uid_usage]
    uid_usage = [d for idx, d in zip(uid_indices, uid_usage) if file.data[idx-1]['str'] == '.'
                                        and file.data[idx+1]['str'] == '(']

    msg = ("Idx unique identifier is not necessarily correlated to array indices.",
           "This function should only be used in DDC's internals.",
           "You probably need code similar to (idx-field_start_idx).value() (.uid())")
    for d in uid_usage:
        report_error(ERROR, file, d['linenr'], msg)

def check_directives(file):
    """
    Test the pragma directives at the start of the file.

    If '#pragma once' is missing from a .hpp file then an error is raised.
    An error is also raised for incorrect include syntax. Correct syntax is when quotes are used to include files
    from the gyselalibxx project and angle brackets are used to include files from other libraries.
    """
    if file.data_xml is None:
        if file.file.suffix != '.hpp':
            return
        try:
            pragma_idx = file.raw.index('pragma')
        except ValueError:
            pragma_idx = None
        if pragma_idx is None or file.raw[pragma_idx-1] != '#' or file.raw[pragma_idx+1] != 'once':
            report_error(FATAL, file, 2, "#pragma once missing from the top of a hpp file")

    if file.file.suffix == '.hpp' and not file.root.findall("./dump/directivelist/directive[@str='#pragma once']"):
        report_error(FATAL, file, 2, "#pragma once missing from the top of a hpp file")

    directives = {d.attrib['linenr']: d.attrib['str'].split() for d in file.root.findall("./dump/directivelist/directive")}
    include_directives = {linenr: words[1] for linenr, words in directives.items() if len(words)>1 and words[0] == '#include'}
    for linenr, include_str in include_directives.items():
        include_file = include_str[1:-1]
        possible_matches = [match for f in global_folders for match in f.glob(f'**/{include_file}')]
        if possible_matches and include_str[0] != '"':
            report_error(STYLE, file, linenr, f'Quotes should be used to include files from the gyselalibxx project ({include_str}->"{include_file}")')
        elif not possible_matches and include_str == '"':
            report_error(STYLE, file, linenr, f'Angle brackets should be used to include files from external libraries ({include_str}-><{include_file}>)')

def check_licence_presence(file):
    """
    Test to ensure that all .hpp files have the license information in the first line.
    """
    if file.file.suffix != '.hpp':
        return
    if file.raw[0]['str'] != '// SPDX-License-Identifier: MIT':
        report_error(ERROR, file, 1, "Licence missing from the top of the file (// SPDX-License-Identifier: MIT)")

if __name__ == '__main__':
    parser = argparse.ArgumentParser("A static analysis scipt to search for common errors using cppcheck")
    parser.add_argument('files', type=str, nargs='*')
    args = parser.parse_args()

    no_file_filter = not args.files
    filter_files = [Path(f).absolute() for f in args.files]

    relevant_files = get_relevant_files()

    multipatch_geom = list(relevant_files.pop('multipatch'))

    cppcheck_command = ['cppcheck', '--dump', '--library=googletest', '--check-level=exhaustive', '--enable=style',
                        '--std=c++17', '--max-ctu-depth=5', '--suppress=unusedStructMember']
    for f in multipatch_geom:
        if no_file_filter or f in filter_files:
            print("------------- Checking ", f, " -------------")
            subprocess.run([*cppcheck_command, f], check=False)

    for geom, files in relevant_files.items():
        if no_file_filter or any(f in filter_files for f in files):
            print("------------- Checking ", geom, " -------------")
            subprocess.run(cppcheck_command + list(files) + [f'--file-filter={f}' for f in filter_files], check=False)

    files = [FileObject(f) for f in HOME_DIR.glob('**/*.dump')]
    for myfile in files:
        check_for_ddc_keywords(myfile)
        search_for_unnecessary_auto(myfile)
        search_for_bad_create_mirror(myfile)
        search_for_bad_aliases(myfile)
        check_struct_names(myfile)
        search_for_uid(myfile)
        search_for_bad_memory(myfile)
        check_directives(myfile)
        check_licence_presence(myfile)

    sys.exit(error_level)
