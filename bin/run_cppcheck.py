#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import subprocess
import shutil
import sys

GYSELALIB_ROOT = Path(__file__).parent.parent

if __name__ == '__main__':
    parser = argparser.ArgumentParser("Use cppcheck to look for common coding errors. This script must be run in a root folder.")
    parser.add_argument('toolchain', 'The toolchain that you would use to compile the code.')
    parser.add_argument('folders', 'The folders in which the code you would like to test is found.')
    parser.add_argument('--disable-style', 'Disable the style checks', action='store_true')
    parser.add_argument('--no-clean-up', 'Disable the clean-up of the generated folder.', action='store_true')
    parser.add_argument('--verbose', 'Run cppcheck in verbose mode', action='store_true')

    args = parser.parse_args()

    toolchain = Path(args.toolchain).absolute()

    build_dir = Path('build_cppcheck').absolute()
    subprocess.run([shutil.which('cmake'), f'-DCMAKE_TOOLCHAIN_FILE=${toolchain}',
               '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON', '-B', build_dir, 'S', '.'],
                   check=True)

    with open(build_dir / 'compile_commands.json', encoding='utf-8') as f:
        config = json.load(f)

    relevant_folders = [r.absolute() for r in args.relevant_folders]

    relevant_configs = [c for c in config if any(r in Path(c['file']).parents for r in relevant_folders)]

    for c in relevant_configs:
        c['command'] = ' '.join(s for s in c['command'].split() if not (s.startswith('-I') and '/vendor/' in s))

    with open(build_dir / 'relevant_compile_commands.json', 'w', encoding='utf-8') as f:
        config = json.dump(relevant_configs, f, indent=2)

    cmd = [shutil.which('cppcheck'), '-i', f'{GYSELALIB_DIR}/vendor',
           '--library=googletest', '--max-ctu-depth=5', '--check-level=exhaustive',
           '--std=c++17', '--suppress=unusedStructMember', '--suppress=useStlAlgorithm',
           '--suppress=knownConditionTrueFalse', '--suppress=ctuOneDefinitionRuleViolation',
           f'--addon={GYSELALIB_DIR}/bin/ci_tools/check_naming_conventions',
           f'--addon={GYSELALIB_DIR}/bin/ci_tools/check_for_memory_misuse',
           f'--project={build_dir}/relevant_compile_commands.json']

    if args.verbose:
        cmd.append('-v')
    if not args.disable_style:
        cmd.append('--enable=style')

    p = subprocess.run(cmd, check=False)

    if not args.no_clean_up:
        shutil.rmtree(build_dir)

    sys.exit(p.exit_code)
