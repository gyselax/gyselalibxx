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
WARNING = "\033[1m\033[33mwarning\033[0m"

possible_error_levels = {WARNING: 0, STYLE: 1, ERROR: 2, FATAL: 3}

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

parallel_functions = ['parallel_for', 'parallel_for_each', 'parallel_transform_reduce']

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
                    for idx in using_indices if self.raw[idx+2]['str'] == '='}

        self.type_names = [{'str': line[1]['str'], 'linenr': line[0]['linenr']} for line in self.lines \
                            if line and line[0]['str'] in ('struct', 'class')]
        for templated_obj in self.root.findall('./dump/TemplateSimplifier/TokenAndName'):
            if Path(templated_obj.attrib['file']) != self._file:
                continue
            linenr = templated_obj.attrib['line']
            self.type_names += [{'str': templated_obj[i+2].attrib['str'] if templated_obj[i+1].attrib['str'] == '...' \
                                        else templated_obj[i+1].attrib['str'],
                                 'linenr': linenr} for i,t in enumerate(templated_obj) \
                                if t.attrib['str'] in ('struct', 'class')]
        self.type_names += [{'str': self.raw[idx+1]['str'], 'linenr': self.raw[idx+1]['linenr']} for idx in using_indices]

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
                       HOME_DIR / 'src' / 'utils' / 'ddc_alias_inline_functions.hpp',
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
    var_data_idx = [file.data.index(v.attrib) for v in auto_variables if v is not None]
    var_names = [v.attrib['str'] for v in auto_variables if v is not None]
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

    if not file.data_xml:
        return

    # Find auto in arguments of lambda functions
    for elem in file.data_xml.findall(".token[@str='[']"):
        start_idx = file.data.index(elem.attrib)+1
        end_idx = next(j for j,a in enumerate(file.data[start_idx:], start_idx) if a['str'] == ']')
        keys = ''.join(a['str'] for a in file.data[start_idx:end_idx]).split(',')
        if not all(len(k) == 0 or k[0] in ('&', '=') for k in keys):
            continue
        end_args_idx = next(j for j,a in enumerate(file.data[end_idx:], end_idx) if a['str'] == ')')
        lambda_args = ' '.join(a['str'] for a in file.data[end_idx+2:end_args_idx]).split(',')
        for a in lambda_args:
            if 'auto' in a:
                report_error(ERROR, file, file.data[start]['linenr'], f"Please use explicit types instead of auto ({a})")

    for elem in chain(file.data_xml.findall(".token[@str='KOKKOS_CLASS_LAMBDA']"), file.data_xml.findall(".token[@str='KOKKOS_LAMBDA']")):
        start_idx = file.data.index(elem.attrib)+3
        end_args_idx = next(j for j,a in enumerate(file.data[start_idx:], start_idx) if a['str'] == ')')
        lambda_args = ' '.join(a['str'] for a in file.data[start_idx:end_args_idx]).split(',')
        for a in lambda_args:
            if 'auto' in a.split():
                report_error(ERROR, file, file.data[start]['linenr'], f"Please use explicit types instead of auto ({a})")

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
        if a_name[0].isupper():
            if val_str.startswith('Idx<') and not a_name.startswith('Idx'):
                report_error(STYLE, file, val[0]['linenr'], f"Index type aliases should start with the keyword Idx ({a_name})")
            elif val_str.startswith('IdxStep<') and not a_name.startswith('IdxStep'):
                report_error(STYLE, file, val[0]['linenr'], f"Index step type aliases should start with the keyword IdxStep ({a_name})")
            elif val_str.startswith('IdxRange<') and not a_name.startswith('IdxRange'):
                report_error(STYLE, file, val[0]['linenr'], f"Index type aliases should start with the keyword IdxRange ({a_name})")

    for alias_key in file.type_names:
        a_name = alias_key['str']
        linenr = alias_key['linenr']
        if any(bad_name in a_name for bad_name in forbidden_substrings):
            report_error(STYLE, file, linenr, f"Name seems to contain DDC keywords. Please use the new naming conventions (using {a_name} = {val_str})")
        if a_name == 'Grid':
            msg = ("The name 'Grid' is too general and should be avoided as it implies that the type "
                   "can be multi-D. Please use Grid1D")
            report_error(STYLE, file, linenr, msg)
        if 'Idx' in a_name and not any(a_name.startswith(valid_start_name) for valid_start_name in ('Idx', 'MultipatchIdx', 'HasIdx', 'InternalIdx')):
            prefix = 'Multipatch' if 'Multipatch' in a_name else ''
            name = a_name.replace('Multipatch','')
            if 'IdxRange' in a_name:
                expected_name = prefix + 'IdxRange' + name.replace('IdxRange','')
            elif 'IdxStep' in a_name:
                expected_name = prefix + 'IdxStep' + name.replace('IdxStep','')
            else:
                expected_name = prefix + 'Idx' + name.replace('Idx','')
            report_error(STYLE, file, linenr, f"Type aliases should start with the type keyword ({a_name} -> {expected_name})")

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
    uid_usage = [d.attrib for d in file.data_xml.findall(".token[@str='uid']") if Path(d.attrib['file']) == file.file]
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

    directives = {d.attrib['linenr']: d.attrib['str'].split() for d in file.root.findall("./dump/directivelist/directive")
                        if Path(d.attrib['file']) == file.file}
    include_directives = {linenr: words[1] for linenr, words in directives.items() if len(words)>1 and words[0] == '#include'}
    for linenr, include_str in include_directives.items():
        include_file = include_str[1:-1]
        possible_matches = [match for f in global_folders for match in f.glob(f'**/{include_file}')]
        if possible_matches and include_str[0] != '"':
            report_error(STYLE, file, linenr, f'Quotes should be used to include files from the gyselalibxx project ({include_str}->"{include_file}")')
        elif not possible_matches and include_str == '"':
            report_error(STYLE, file, linenr, f'Angle brackets should be used to include files from external libraries ({include_str}-><{include_file}>)')

def update_aliases(all_files):
    """
    Update the aliases using information from geometry.hpp and species_info.hpp. This ensures that as much information
    as possible is available to describe the type. E.g. this helps identify the memory space of a Field.
    """
    if len(all_files) == 0:
        return

    spec_file = next(f for f in all_files if f.file == HOME_DIR / 'src/speciesinfo/species_info.hpp')
    geom_files = [f for f in all_files if f.file.name == 'geometry.hpp' or 'geometries' in f.file.parts]
    multipatch_geom_file_names = {f.file.name: f for f in geom_files if f.file.name != 'geometry.hpp'}

    for file in all_files:
        if file.data_xml is None:
            continue
        directives = {d.attrib['linenr']: d.attrib['str'].split() for d in file.root.findall("./dump/directivelist/directive")
                            if Path(d.attrib['file']) == file.file}
        include_directives = {linenr: words[1] for linenr, words in directives.items() if len(words)>1 and words[0] == '#include'}
        for linenr, include_str in include_directives.items():
            include_file = include_str[1:-1]
            if include_file == 'species_info.hpp':
                file.aliases.update(spec_file.aliases)
            elif include_file == 'geometry.hpp':
                relevant_geom = next(p for p in file.file.parts if p.startswith('geometry'))
                geom_file = next(f for f in geom_files if relevant_geom in f.file.parts)
                file.aliases.update(geom_file.aliases)
            elif include_file in multipatch_geom_file_names:
                file.aliases.update(multipatch_geom_file_names[include_file].aliases)


def check_kokkos_lambda_use(file):
    """
    Check that KOKKOS_LAMBDA expressions are present and that they are correctly used.

    The checks are:
    - Check that lambda functions passed to ddc::parallel_X functions use KOKKOS_LAMBDA or KOKKOS_CLASS_LAMBDA
    - Check that class variables are not used in lambda functions passed to ddc::parallel_X functions.
    """
    if not file.data_xml:
        return

    for p in parallel_functions:
        for elem in file.data_xml.findall(f".token[@str='{p}']"):
            idx = file.data.index(elem.attrib)
            open_brackets = file.data[idx+1]
            close_brackets_id = open_brackets['link']
            scope = open_brackets['scope']
            end_idx = file.data.index(file.data_xml.find(f".token[@id='{close_brackets_id}'][@scope='{scope}']").attrib)
            arg_splitters = [idx+1]
            i = idx+3
            while i < end_idx:
                if file.data[i]['str'] == ',':
                    arg_splitters.append(i)
                    i+=1
                elif file.data[i]['str'] in ('(','{', '['):
                    close_brackets_id = file.data[i]['link']
                    scope = file.data[i]['scope']
                    i = file.data.index(file.data_xml.find(f".token[@id='{close_brackets_id}'][@scope='{scope}']").attrib)+1
                else:
                    i+=1

            if file.data[arg_splitters[-1]+1]['str'] == '[':
                last_arg = ' '.join(a['str'] for a in file.data[arg_splitters[-1]+1:end_idx])
                correct_last_arg = 'KOKKOS_LAMBDA '+last_arg[last_arg.find(']')+1:]
                report_error(ERROR, file, elem.attrib['linenr'], f'The lambda function passed to the function ddc::{p} must be a KOKKOS_LAMBDA\n{last_arg}\nShould be:\n{correct_last_arg}')

            lambda_body_start = next(i for i,a in enumerate(file.data[arg_splitters[-1]+1:], arg_splitters[-1]+1) if a['str'] == '{')
            lambda_body_end_id = file.data[lambda_body_start]['link']
            lambda_body_scope = file.data[lambda_body_start]['scope']
            end_body_idx = file.data.index(file.data_xml.find(f".token[@id='{lambda_body_end_id}'][@scope='{lambda_body_scope}']").attrib)
            for a in file.data[lambda_body_start+1:end_body_idx]:
                var_id = a.get('variable', None)
                if var_id:
                    var = next(v for v in file.variables if v['id'] == var_id)
                    var_scope = var['scope']
                    scope_info = file.root.find(f"./dump/scopes/scope[@id='{var_scope}']").attrib
                    if scope_info['type'] == 'Class':
                        report_error(ERROR, file, a['linenr'], f'Please create a local variable with the suffix "_proxy" to store the class variable that is used in a parallel loop ({a["str"]})')

def check_licence_presence(file):
    """
    Test to ensure that all .hpp files have the license information in the first line.
    """
    if file.file.suffix != '.hpp':
        return
    if file.raw[0]['str'] != '// SPDX-License-Identifier: MIT':
        report_error(ERROR, file, 1, "Licence missing from the top of the file (// SPDX-License-Identifier: MIT)")

def check_exec_space_usage(file):
    """
    Test to ensure that the operator() of a Field or FieldMem is not used from the wrong space (CPU/GPU)
    """
    # Treat each configuration separately
    for config in file.root.findall('./dump'):
        tokenlist = [t.attrib for t in config.find('./tokenlist')]

        # Find all scopes and specify DefaultHostExecutionSpace by default
        scopes = {s.attrib['id'] : s.attrib for s in config.findall(".//scope")}
        for s in scopes.values():
            s['exec_space'] = 'DefaultHostExecutionSpace'

        known_execution_spaces = ('DefaultHostExecutionSpace', 'DefaultExecutionSpace')

        # Find all markers indicating a Kokkos scope and note the correct execution space
        kokkos_scope_markers = chain(config.findall(".//token[@str='KOKKOS_LAMBDA']"),
                                     config.findall(".//token[@str='KOKKOS_FUNCTION']"),
                                     config.findall(".//token[@str='KOKKOS_INLINE_FUNCTION']"))

        for m in kokkos_scope_markers:
            # Identify scope id
            current_idx = tokenlist.index(m.attrib)
            start_idx = next(i for i,t in enumerate(tokenlist[current_idx:],current_idx) if t['str'] =='{' or t['str'] == ';')
            if tokenlist[start_idx]['str'] == ';':
                continue

            scope_id = tokenlist[start_idx+1]['scope']
            scope = scopes[scope_id]

            if tokenlist[current_idx]['str'] == 'KOKKOS_LAMBDA':
                ast_parent = config.find(f".//token[@id='{tokenlist[current_idx-1]['astParent']}']").attrib
                parallel_func_start_idx = tokenlist.index(ast_parent)
                args = [t['str'] for t in tokenlist[parallel_func_start_idx+1:current_idx]]
                if 'DefaultHostExecutionSpace' in args:
                    continue

            # Update execution space
            scopes[scope_id]['exec_space'] = 'DefaultExecutionSpace'

            # Ensure that any nested scopes are also updated
            to_pass_down = set([scope_id])
            while to_pass_down:
                new_to_pass_down = set()
                for s in scopes.values():
                    if s.get('nestedIn', None) in to_pass_down:
                        new_to_pass_down.add(s['id'])
                        s['exec_space'] = 'DefaultExecutionSpace'
                to_pass_down = new_to_pass_down

        field_types = ('Field', 'FieldMem', 'DField', 'DFieldMem', 'ConstField', 'DConstField')

        # Find all variables that may be defined on either CPU or GPU
        field_variables = {}
        auto_variables = {}
        for v in config.findall("./variables/var"):
            # Identify the type
            type_start_token = v.get('typeStartToken')
            type_end_token = v.get('typeEndToken')
            name_token = v.get('nameToken')
            type_start_idx = tokenlist.index(config.find(f"./tokenlist/token[@id='{type_start_token}']").attrib)
            type_end_idx = tokenlist.index(config.find(f"./tokenlist/token[@id='{type_end_token}']").attrib)
            type_elements = [t['str'] for t in tokenlist[type_start_idx:type_end_idx+1]]
            name = config.find(f"./tokenlist/token[@id='{name_token}']").attrib['str'] if name_token != '0' else ''
            type_descr = ' '.join(type_elements)
            if type_descr in file.aliases:
                type_elements = [t['str'] for t in file.aliases[type_descr]]
                type_descr = ' '.join(type_elements)

            # If the type is a Field or a FieldMem then identify the space where the memory is located
            if (type_elements[0] in field_types) or (type_elements[0] == 'host_t' and type_elements[1] in field_types):
                attribs = v.attrib
                attribs['type'] = type_descr
                field_variables[name] = attribs
                if 'host_t' == type_elements[0]:
                    attribs['exec_space'] = 'DefaultHostExecutionSpace'
                elif 'DefaultExecutionSpace' in type_elements:
                    attribs['exec_space'] = 'DefaultExecutionSpace'
                elif 'DefaultHostExecutionSpace' in type_elements:
                    attribs['exec_space'] = 'DefaultHostExecutionSpace'
                else:
                    if '<' in type_descr:
                        start_template_args = type_descr.index('<')
                        end_template_args = type_descr.rindex('>')
                        template_args = type_descr[start_template_args+1:end_template_args].split(' ')
                        args = []
                        last_idx = 0
                        template_depth = 0
                        for i,a in enumerate(template_args):
                            if a == '<':
                                template_depth += 1
                            elif a == '>':
                                template_depth -= 1
                            elif template_depth == 0 and a == ',':
                                args.append(' '.join(template_args[last_idx:i]))
                                last_idx = i+1
                        args.append(' '.join(template_args[last_idx:]))
                        if type_elements[0].startswith('D'):
                            args.insert(0,'double')
                        if 'Mem' in type_descr and len(args) > 2:
                            if 'DeviceAllocator' in args[-1]:
                                attribs['exec_space'] = 'DefaultExecutionSpace'
                            elif 'HostAllocator' in args[-1]:
                                attribs['exec_space'] = 'DefaultHostExecutionSpace'
                            else:
                                attribs['exec_space'] = args[-1].split(',')[-1].strip(' >')
                        elif 'Mem' not in type_descr and len(args) > 3:
                            attribs['exec_space'] = args[-1]
                        else:
                            attribs['exec_space'] = 'DefaultExecutionSpace'
                    else:
                        attribs['exec_space'] = 'DefaultExecutionSpace?'
            # If a type is auto then further checks are required
            elif type_descr == 'auto':
                attribs = v.attrib
                attribs['type'] = type_descr
                auto_variables[name] = attribs

        # Check auto variables to identify Field/FieldMem types
        for name, v in auto_variables.items():
            v_id = v['id']
            var_usage = config.findall(f"./tokenlist/token[@variable='{v_id}']")
            for u in var_usage:
                idx = tokenlist.index(u.attrib)
                if tokenlist[idx+1]['str'] == '=':
                    phrase_end = next(i for i,t in enumerate(tokenlist[idx+2:],idx+2) if t['str'] == ';')
                    keys = [t['str'] for t in tokenlist[idx+2:phrase_end]]
                    # If the auto variable is created by a call to a mirror function then identify the memory space
                    if any(k in mirror_functions for k in keys):
                        if 'DefaultExecutionSpace' in keys:
                            v['exec_space'] = 'DefaultExecutionSpace'
                        elif 'DefaultHostExecutionSpace' in keys:
                            v['exec_space'] = 'DefaultHostExecutionSpace'
                        else:
                            v['exec_space'] = 'DefaultHostExecutionSpace'
                        field_variables[name] = v

        # For each Field/FieldMem type check that it's contents are only accessed on the correct memory space
        for name, v in field_variables.items():
            v_id = v['id']
            var_usage = config.findall(f"./tokenlist/token[@variable='{v_id}']")
            for u in var_usage:
                idx = tokenlist.index(u.attrib)
                scope_id = tokenlist[idx]['scope']
                scope = scopes[scope_id]
                if tokenlist[idx]['file'] != str(file.file):
                    continue
                # Exclude initialisation where the variable is followed by ( but the operator() is not called
                if tokenlist[idx-1]['id'] == v['typeEndToken'] or scope['type'] in ('Class', 'Struct'):
                    continue
                # If accessing the contents then we must be on the correct execution space
                if tokenlist[idx+1]['str'] == '(':
                    exec_space = scope['exec_space']
                    memory_space = v['exec_space']
                    level = None

                    if exec_space in known_execution_spaces and memory_space in known_execution_spaces and exec_space != memory_space:
                        level = FATAL
                    elif exec_space == 'DefaultHostExecutionSpace' and exec_space != memory_space:
                        level = WARNING
                        if file.file.name == 'species_info.hpp':
                            # Deactivate warning for species_info constructor which purposely has different constructors for different spaces
                            level = None

                    if level:
                        msg = ("Attempted memory access from the wrong execution space. "
                               f"The current execution space is \"{scope['exec_space']}\" but the variable "
                               f"{name} is stored on a space compatible with \"{v['exec_space']}\"")
                        report_error(level, file, tokenlist[idx]['linenr'], msg)

                    if scope['exec_space'] == 'DefaultExecutionSpace' and (v['type'] == 'auto' or 'FieldMem' in v['type']):
                        msg = ''
                        if v['type'] == 'auto':
                            msg = "An auto type may be either a Field or a FieldMem depending on compilation paramaters. "
                        msg += ("The FieldMem copy operator is deleted to avoid accidental memory allocation. "
                                "You probably do not want a copy of the FieldMem allocated on each GPU thread. "
                                "Please use a Field to access the data that has already been allocated.")
                        report_error(FATAL, file, tokenlist[idx]['linenr'], msg)

        for s_id, scope in scopes.items():
            if scope['exec_space'] == 'DefaultHostExecutionSpace':
                continue
            relevant_code = config.findall(f".//token[@scope='{s_id}']")
            code_keys = [c.attrib['str'] for c in relevant_code]
            if 'std' in code_keys:
                idx = code_keys.index('std')
                if code_keys[idx+1] == '::':
                    msg = "Std functions are not designed to run on GPU. You may wish to check if there is an equivalent Kokkos:: function?"
                    report_error(FATAL, file, relevant_code[idx].attrib['linenr'], msg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("A static analysis scipt to search for common errors using cppcheck")
    parser.add_argument('files', type=str, nargs='*')
    args = parser.parse_args()

    no_file_filter = not args.files
    filter_files = [Path(f).absolute() for f in args.files]
    spec_info = HOME_DIR / 'src/speciesinfo/species_info.hpp'

    relevant_files = get_relevant_files()

    multipatch_geom = list(relevant_files.pop('multipatch'))

    cppcheck_command = ['cppcheck', '--dump', '--library=googletest', '--check-level=exhaustive', '--enable=style',
                        '--std=c++17', '--max-ctu-depth=5', '--suppress=unusedStructMember', '--suppress=useStlAlgorithm',
                        '--error-exitcode=1']
    for f in multipatch_geom:
        if no_file_filter or f in (*filter_files, spec_info):
            print("------------- Checking ", f, " -------------")
            p = subprocess.run([*cppcheck_command, f], check=False)
            if p.returncode:
                error_level = max(error_level, possible_error_levels[STYLE])

    for geom, files in relevant_files.items():
        if no_file_filter or any(f in filter_files for f in files):
            print("------------- Checking ", geom, " -------------")
            p = subprocess.run(cppcheck_command + list(files) + [f'--file-filter={f}' for f in filter_files+['*geometry.hpp', str(spec_info)]], check=False)
            if p.returncode:
                error_level = max(error_level, possible_error_levels[STYLE])

    files = [FileObject(f) for f in HOME_DIR.glob('**/*.dump')]
    for myfile in files:
        check_for_ddc_keywords(myfile)
        search_for_unnecessary_auto(myfile)
        search_for_bad_create_mirror(myfile)
        search_for_bad_aliases(myfile)
        check_struct_names(myfile)
        search_for_uid(myfile)
        check_directives(myfile)
        check_licence_presence(myfile)
        check_kokkos_lambda_use(myfile)

    update_aliases(files)

    for myfile in files:
        search_for_bad_memory(myfile)
        check_exec_space_usage(myfile)

    sys.exit(error_level)
