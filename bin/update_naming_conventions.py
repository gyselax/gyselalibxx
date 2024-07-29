import argparse
from itertools import chain
from pathlib import Path
import re
import shutil
import subprocess

HOME_DIR = Path(__file__).parent.parent

def update_names(old_name, new_name, applicable_files, complete_word = True):
    """
    Update the name of an object across the files.

    Parameters
    ----------
    old_name : str
        The current name found in the code.
    new_name : str
        The name which should be used instead.
    applicable_files : list[Path]
        A list of the files whose contents should be modified.
    complete_word : bool, default=True
        True if only complete words should be replaced, False for complex regex expressions.
    """
    if complete_word:
        old_name = rf'\<{old_name}\>'
    subprocess.run(['sed', '-i', f's/{old_name}/{new_name}/g'] + applicable_files)

key_map = {'<':'>',
           '(':')',
           '[':']',
           '{':'}'}

def find_argument(open_tag, expr):
    """
    Find the arguments at the start of an expression.

    The arguments may be template arguments, function arguments etc.
    The function takes a string beginning after the opening tag and
    looks for the closing tag which matches the opening tag. It then
    returns the arguments between the two tags.
    """
    close_tag = key_map[open_tag]
    opened = [None]
    i = 0
    while len(opened):
        if expr[i] == open_tag:
            opened.append(i)
        elif expr[i] == close_tag:
            opened.pop()

        i+= 1
    return expr[:i]

def get_new_Coord_name(alias_match):
    """
    Get the new name of a coordinate from a regex describing a type alias (see remove_conflicting_names).
    """
    groups = alias_match.groups()
    if 'continuous_element_type' in groups[0]:
        return 'Coord1D'
    else:
        tags = groups[2]
        if '...' in tags:
            return 'CoordND'
        else:
            ntags = tags.count(',')+1
            return f'Coord{ntags}D'

def remove_conflicting_names(applicable_files):
    """
    Remove names used in the code which conflict with the new gysela templated types.
    The names are:
    -  Coord

    Parameters
    ----------
    applicable_files : list[Path]
        A list of the files whose contents should be modified.
    """
    for file in applicable_files:
        with open(file, 'r') as f:
            code = f.read()

        type_aliases = list(re.finditer('using Coord = ((ddc::)?Coordinate<([a-zA-Z0-9_:, \.]+)>|[a-zA-Z0-9 ]+::continuous_element_type)', code))
        if not type_aliases:
            continue

        if len(type_aliases) == 1:
            new_name = get_new_Coord_name(type_aliases[0])
            code = re.sub(r'\bCoord\b(?!<)', new_name, code)
        else:
            alias_match = type_aliases[0]
            while alias_match:
                scope_start = code.rfind('{', 0, alias_match.start())
                scope_contents = find_argument('{', code[scope_start + 1:])
                scope_end = len(scope_contents) + scope_start + 1

                new_name = get_new_Coord_name(alias_match)
                code = code[:scope_start + 1] + re.sub(r'\bCoord\b(?!<)', new_name, scope_contents) \
                        + code[scope_end:]

                alias_match = re.search('using Coord = ((ddc::)?Coordinate<([a-zA-Z0-9_:, \.]+)>|[a-zA-Z0-9 ]+::continuous_element_type)', code)

        with open(file, 'w') as f:
            f.write(code)
            if not code.endswith('\n'):
                f.write('\n')

def replace_array(old_name, new_name, code, memory_arg_position):
    """
    Replace the name of an array-like object, inserting host_t if necessary.

    Parameters
    ----------
    old_name : str
        The current name found in the code.
    new_name : str
        The name which should be used instead.
    code : str
        The code to be modified.
    memory_arg_position : int
        The index indicating the position where the memory space is specified
        in the template arguments.
    """
    location = re.search(f'{old_name}<', code)
    while location:
        arg = find_argument('<', code[location.end():])
        if arg.count(',') < memory_arg_position:
            code = code[:location.start()] + f'host_t<{new_name}<' + arg + '>' + code[location.end() + len(arg):]
        else:
            code = code[:location.start()] + f'{new_name}<' + arg + code[location.end() + len(arg):]
        location = re.search(f'{old_name}<', code)

    return code

def replace_device_array(old_name, new_name, code):
    """
    Replace the name of an array-like object found on the device.
    The device tag is removed as it is specified by default in Gysela naming conventions.

    Parameters
    ----------
    old_name : str
        The current name found in the code.
    new_name : str
        The name which should be used instead.
    code : str
        The code to be modified.
    memory_arg_position : int
        The index indicating the position where the memory space is specified
        in the template arguments.
    """
    location = re.search(rf'device_t<\s*{old_name}<', code)
    while location:
        arg = find_argument('<', code[location.end():])
        close_device = code.find('>', location.end() + len(arg))
        code = code[:location.start()] + f'{new_name}<' + arg + code[close_device+1:]
        location = re.search(f'device_t<{old_name}<', code)

    return code

def replace_memory_block(old_name, new_name, code, memory_arg_position):
    """
    Replace the name of an array-like object, inserting and removing host_t/device_t if necessary.

    Parameters
    ----------
    old_name : str
        The current name found in the code.
    new_name : str
        The name which should be used instead.
    code : str
        The code to be modified.
    memory_arg_position : int
        The index indicating the position where the memory space is specified
        in the template arguments.
    """
    code = replace_device_array(old_name, new_name, code)
    code = replace_array(old_name, new_name, code, memory_arg_position)
    return re.sub(f'\<{old_name}\>', new_name, code)

def replace_domain_variable_name(code):
    """
    Replace the word domain or dom in the code. Only places where the word is
    proceeded and succeeded by either an underscore or a space are modified.

    Parameters
    ----------
    code : str
        The code to be modified.
    """
    code = re.sub('(_|\b)(dom)(_|\b)', r'\1idx_range\3', code)

    location = re.search('domain', code)
    while location:
        line_start_idx = code.rfind('\n', 0, location.start())
        last_comment_block_start = code.rfind('/*', 0, location.start())
        last_comment_block_end = code.rfind('*/', 0, location.start())
        in_comment = '//' in code[line_start_idx:location.start()] or \
                        last_comment_block_start > last_comment_block_end

        new_text = 'index range' if in_comment else 'idx_range'
        code = code[:location.start()] + new_text + code[location.end():]

        location = re.search(rf'domain', code)

    # Revert for bad renames
    code = code.replace('@param index range', '@param idx_range')
    code = code.replace('@param[in] index range', '@param[in] idx_range')

    # Revert for DDC names
    code = code.replace('::get_idx_range', '::get_domain')
    code = code.replace('discrete_idx_range_type', 'discrete_domain_type')
    code = code.replace('template get_idx_range', 'template get_domain')
    code = code.replace('batched_spline_idx_range', 'batched_spline_domain')

    return code.replace(' a index ', ' an index ')

def update_names_using_context_clues(applicable_files):
    """
    Update the names of objects when the new name is dependant on context clues.

    This must be done via parsing to find those clues.
    The function updates the names of Chunks, Spans and Views. These require
    context clues to handle memory spaces correctly.
    It also updates variable names describing domains
    """
    for fname in applicable_files:
        with open(fname, 'r') as f:
            orig_code = f.read()

        code = replace_memory_block('ddc::Chunk', 'FieldMem', orig_code, 2)
        code = replace_memory_block('ddc::ChunkSpan', 'Field', code, 3)
        code = replace_memory_block('ddc::ChunkView', 'ConstField', code, 3)

        code = replace_domain_variable_name(code)

        if code != orig_code:
            with open(fname, 'w') as f:
                f.write(code)
                if not code.endswith('\n'):
                    f.write('\n')

def get_geometry_name_shortcuts(geometry_file):
    """
    Read the geometry.hpp file and look for all types/aliases that should be renamed to
    follow the new conventions.

    Parameters
    ----------
    geometry_file : str
        The path to the geometry file.

    Returns
    -------
    dict[str, str]
        A dictionary whose keys are types to be renamed and whose values are the new
        names.
    """
    with open(geometry_file, 'r') as f:
        code = f.read()

    structures = re.finditer('(struct|using|class) ([a-zA-Z0-9]+)', code)

    name_map = {}

    species_names = {'IdxSp', 'IdxRangeSp', 'IdxStepSp', 'FieldMemSp', 'DFieldMemSp', 'IFieldSp', 'ConstFieldSp', 'DConstFieldSp', 'FieldSp', 'DFieldSp'}

    for s in structures:
        struct_name = s[2]
        if struct_name in species_names:
            continue
        if struct_name == 'RDim':
            name_map[struct_name] = 'Dim'
        elif struct_name in ('DDim', 'IDim'):
            name_map[struct_name] = 'Grid1D'
        elif struct_name.startswith('RDim'):
            new_name = struct_name[4:]
            if not new_name[0].isalpha():
                new_name = 'Dim' + new_name
            name_map[struct_name] = new_name
        elif struct_name.startswith('IDim') or struct_name.startswith('DDim'):
            name_map[struct_name] = 'Grid' + struct_name[4:]
        elif struct_name.startswith('Index'):
            name_map[struct_name] = 'Idx' + struct_name[5:]
        elif struct_name.startswith('IVect'):
            name_map[struct_name] = 'IdxStep' + struct_name[5:]
        elif struct_name.startswith('DVec'):
            name_map[struct_name] = 'IdxStep' + struct_name[4:]
        elif struct_name.startswith('IDomain'):
            name_map[struct_name] = 'IdxRange' + struct_name[7:]
        elif struct_name.startswith('BSDomain'):
            name_map[struct_name] = 'BSIdxRange' + struct_name[8:]
        elif struct_name.endswith('Index'):
            name_map[struct_name] = struct_name[:-5] + 'Idx'
        elif struct_name.startswith('Field'):
            name_map[struct_name] = 'FieldMem' + struct_name[5:]
        elif struct_name.startswith('Span'):
            name_map[struct_name] = 'Field' + struct_name[4:]
        elif struct_name.startswith('View'):
            name_map[struct_name] = 'ConstField' + struct_name[4:]
        elif struct_name.startswith('DField'):
            name_map[struct_name] = 'DFieldMem' + struct_name[6:]
        elif struct_name.startswith('DSpan'):
            name_map[struct_name] = 'DField' + struct_name[5:]
        elif struct_name.startswith('DView'):
            name_map[struct_name] = 'DConstField' + struct_name[5:]
        elif 'DomainType' in struct_name:
            name_map[struct_name] = struct_name.replace('DomainType', 'IdxRange')
        elif 'Domain' in struct_name:
            name_map[struct_name] = struct_name.replace('Domain', 'IdxRange')
        elif 'DDom' in struct_name:
            name_map[struct_name] = struct_name.replace('DDom', 'IdxRange')
        elif 'ChunkSpan' in struct_name:
            name_map[struct_name] = struct_name.replace('ChunkSpan', 'Field')
        elif 'ChunkView' in struct_name:
            name_map[struct_name] = struct_name.replace('ChunkView', 'ConstField')
        elif 'Span' in struct_name:
            name_map[struct_name] = struct_name.replace('Span', 'Field')
        elif 'View' in struct_name:
            name_map[struct_name] = struct_name.replace('View', 'ConstField')
        elif struct_name.startswith('DElem'):
            name_map[struct_name] = 'Idx' + struct_name[5:]

    return name_map

def add_missing_include():
    """
    Add the include to get access to the new aliases in files that have been modified to use these aliases.
    Additionally add the necessary linkage in the CMakeLists.txt file.
    """
    p = subprocess.run(['git', 'diff', '--name-only'], capture_output=True, cwd=HOME_DIR, universal_newlines=True)
    modified_files = [HOME_DIR / file for file in p.stdout.split() if file[-4:] in ('.hpp', '.cpp')]
    for file in modified_files:
        subprocess.run(['sed', '-i', '-e', "/#pragma once/a\\", '-e', '#include "ddc_aliases.hpp"', file])

    folders = set(file.parent for file in modified_files)
    try:
        folders.remove(Path(HOME_DIR / 'src' / 'utils'))
    except KeyError:
        pass

    for f in folders:
        cmake = f / 'CMakeLists.txt'
        with open(cmake, 'r') as f:
            code = f.read()

        key = 'target_link_libraries('
        key_len = len(key)
        lib_loc = code.find(key)
        while lib_loc != -1:
            lib_loc += key_len
            args = find_argument('(', code[lib_loc:])
            if 'gslx::utils' not in args:
                last_gslx_lib = args.rfind('gslx::')
                if last_gslx_lib == -1:
                    last_gslx_lib = args.rfind(' ')-1
                line_start = args.rfind('\n', 0, last_gslx_lib)
                line_end = args.find('\n', last_gslx_lib)
                code = code[:lib_loc] + args[:line_end] + args[line_start:last_gslx_lib] + \
                        'gslx::utils\n' + args[line_end:] + code[lib_loc+len(args):]
            lib_loc = code.find(key, lib_loc)

        with open(cmake, 'w') as f:
            f.write(code)
            if not code.endswith('\n'):
                f.write('\n')

def introduce_get_idx_range(applicable_files):
    """
    Replace calls to the DDC .domain() method with calls to get_idx_range.

    Parameters
    ----------
    applicable_files : list[Path]
        A list of the files whose contents should be modified.
    """
    update_names(r'\([a-zA-Z0-9_.]\+\)\.spline_domain()', r'get_spline_idx_range(\1)', applicable_files, complete_word = False)
    update_names(r'\([a-zA-Z0-9_.]\+\)\.domain()', r'get_idx_range(\1)', applicable_files, complete_word = False)
    update_names(r'\([a-zA-Z0-9_.]\+\)\.domain<\([a-zA-Z0-9, .]*\)>()', r'get_idx_range<\2>(\1)', applicable_files, complete_word = False)
    update_names(r'\([a-zA-Z0-9_.]\+\)\.get_domain()', r'get_idx_range(\1)', applicable_files, complete_word = False)
    update_names(r'\([a-zA-Z0-9_.]\+\)\.get_domain<\([a-zA-Z0-9, .]*\)>()', r'get_idx_range<\2>(\1)', applicable_files, complete_word = False)
    update_names(r'ddc::get_domain(\([a-zA-Z0-9_]\+\))', r'get_idx_range(\1)', applicable_files, complete_word = False)
    update_names(r'ddc::get_domain<\([a-zA-Z0-9, .]*\)>(\([a-zA-Z0-9_]\+\))', r'get_idx_range<\1>(\2)', applicable_files, complete_word = False)
    update_names(r'ddc::domain(\([a-zA-Z0-9_]\+\))', r'get_idx_range(\1)', applicable_files, complete_word = False)
    update_names(r'ddc::domain<\([a-zA-Z0-9, .]*\)>(\([a-zA-Z0-9_]\+\))', r'get_idx_range<\1>(\2)', applicable_files, complete_word = False)

def introduce_field_getters(applicable_files):
    """
    Replace calls to the DDC .domain() method with calls to get_idx_range.

    Parameters
    ----------
    applicable_files : list[Path]
        A list of the files whose contents should be modified.
    """
    update_names(r'\([a-zA-Z0-9_.]\+\)\.span_view()', r'get_field(\1)', applicable_files, complete_word = False)
    update_names(r'\([a-zA-Z0-9_.]\+\)\.span_cview()', r'get_const_field(\1)', applicable_files, complete_word = False)

def update_local_type_aliases(applicable_hpp_files, applicable_cpp_files):
    """
    Loop over files to find type aliases which are local to classes or mains.
    """
    treated_cpp_files = set()
    for file in applicable_hpp_files:
        name_map = get_geometry_name_shortcuts(file)
        files_to_update = [file]
        cpp_file = file.with_suffix('.cpp')
        if cpp_file.exists():
            files_to_update.append(cpp_file)
            treated_cpp_files.add(cpp_file)

        for old_name, new_name in name_map.items():
            update_names(old_name, new_name, files_to_update)

    for file in set(applicable_cpp_files).difference(treated_cpp_files):
        name_map = get_geometry_name_shortcuts(file)

        for old_name, new_name in name_map.items():
            update_names(old_name, new_name, [file])

if __name__ == '__main__':
    parser = argparse.ArgumentParser("A script to automatise some of the renaming. This script is not perfect so please check the changes and the compilation")
    parser.add_argument('folders', type=Path, help='The folder where files should be renamed', nargs='+')
    args = parser.parse_args()

    folders = [f.resolve().absolute() for f in args.folders]
    for folder in folders:
        assert HOME_DIR in folder.parents

    src_folder = (HOME_DIR / 'src')
    test_folder = (HOME_DIR / 'tests')
    simu_folder = (HOME_DIR / 'simulations')

    # Exclude the file which renames DDC aliases
    hpp_files = [h for folder in folders for h in folder.glob(f'**/*.hpp') \
                    if h.stem not in ('ddc_aliases', 'ddc_helper')]
    main_files = []
    for folder in folders:
        # Files which may contain a main
        if folder in (test_folder, simu_folder) or test_folder in folder.parents or simu_folder in folder.parents:
            main_files.extend(folder.glob('**/*.cpp'))
    cpp_files = [h for folder in folders for h in folder.glob('**/*.cpp')]
    files_to_update = hpp_files + cpp_files

    remove_conflicting_names(chain(hpp_files, main_files))

    # Update basic names with no additional considerations
    update_names('DiscreteSubDomain', 'IdxRangeSlice', files_to_update)
    update_names('ddc::Coordinate', 'Coord', files_to_update)
    update_names('ddc::DiscreteElement', 'Idx', files_to_update)
    update_names('ddc::DiscreteVector', 'IdxStep', files_to_update)
    update_names('ddc::DiscreteDomain', 'IdxRange', files_to_update)
    update_names('ddc::NonUniformPointSampling', 'NonUniformGridBase', files_to_update)
    update_names('ddc::UniformPointSampling', 'UniformGridBase', files_to_update)
    update_names('IDimSp', 'Species', files_to_update)
    update_names('discrete dimension', 'discretised dimension', files_to_update)

    introduce_get_idx_range(files_to_update)
    introduce_field_getters(files_to_update)

    # Update chunk names taking memory space into consideration
    update_names_using_context_clues(files_to_update)

    add_missing_include()

    name_map = get_geometry_name_shortcuts(HOME_DIR / 'src/speciesinfo/species_info.hpp')

    # Update species names
    for old_name, new_name in name_map.items():
        update_names(old_name, new_name, files_to_update)

    for folder in folders:
        for geom_file in folder.glob('**/geometry.hpp'):
            name_map = get_geometry_name_shortcuts(geom_file)

            geom = geom_file.parents[1]

            # Get files associated with the geometry
            files_to_update = list(chain(geom.glob('**/*.hpp'), geom.glob('**/*.cpp'),
                                    (HOME_DIR / 'tests' / geom.stem).glob('**/*.hpp'),
                                    (HOME_DIR / 'tests' / geom.stem).glob('**/*.cpp'),
                                    (HOME_DIR / 'simulations' / geom.stem).glob('**/*.hpp'),
                                    (HOME_DIR / 'simulations' / geom.stem).glob('**/*.cpp'),
                                    ))

            # Update names within the geometry
            for old_name, new_name in name_map.items():
                update_names(old_name, new_name, files_to_update)

    update_local_type_aliases([f for f in hpp_files if not set(f.parts).intersection({'data_types', 'geometry.hpp', 'species_info.hpp'})], cpp_files)

    clang_format = shutil.which('clang-format')
    if not clang_format:
        clang_format = Path('/data/gyselarunner/clang-format')
        if not clang_format.exists():
            print("WARNING : Clang-format not found. Exiting without cleaning indenting")

    for f in chain(hpp_files, cpp_files):
        subprocess.run([clang_format, '-i', f])
