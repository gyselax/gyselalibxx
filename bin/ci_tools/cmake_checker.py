"""
A script to ensure that any rules about CMake contents are respected.
"""
import glob
from pathlib import Path
import sys

HOME_DIR = Path(__file__).parent.parent.parent.absolute()

if __name__ == '__main__':
    cmake_files = glob.glob(str(HOME_DIR / '**' / 'CMakeLists.txt'), recursive=True)

    success = True

    for f_name in cmake_files:
        with open(f_name, encoding='utf-8') as f:
            code = f.read()

        # Check that all tests are discovered in PRE_TEST mode
        while 'gtest_discover_tests' in code:
            code = code.split('gtest_discover_tests', 1)[1].split('(', 1)[1]
            gtest_arg_str, code = code.split(')', 1)
            gtest_args = gtest_arg_str.split()
            try:
                discovery_index = gtest_args.index('DISCOVERY_MODE')
            except ValueError:
                discovery_index = None
            if discovery_index is None or gtest_args[discovery_index+1] != 'PRE_TEST':
                success = False
                print(f"Tests for target {gtest_args[0]} should be discovered before running the test (not during build). Please add 'DISCOVERY_MODE PRE_TEST' to the arguments of gtest_discover_tests in {f_name}")

    if not success:
        sys.exit(1)
