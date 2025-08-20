"""
Check that the naming conventions documented in CODING_STANDARD.md are respected.
This file should be run using ./bin/run_cppcheck
"""
from pathlib import Path
import re
import sys

import cppcheckdata #pylint: disable=import-error

ddc_keyword_map = {'DiscreteElement': 'Idx',
                   'DiscreteVector' : 'IdxStep',
                   'DiscreteDomain' : 'IdxRange',
                   'Chunk'          : 'FieldMem',
                   'ChunkSpan'      : 'Field',
                   'ChunkView'      : 'ConstField',
                   'Coordinate'     : 'Coord',
                   'span_view'      : 'get_field',
                   'span_cview'     : 'get_const_field',
                   'domain'         : 'get_idx_range',
                   'spline_domain'  : 'get_spline_idx_range',
                   'UniformPointSampling' : 'UniformGridBase',
                   'NonUniformPointSampling' : 'NonUniformGridBase'}

HOME_DIR = Path(__file__).parent.parent.parent.absolute()

file_exceptions = [HOME_DIR / 'src' / 'utils' / 'ddc_helper.hpp',
                   HOME_DIR / 'src' / 'utils' / 'ddc_aliases.hpp',
                   HOME_DIR / 'src' / 'utils' / 'ddc_alias_inline_functions.hpp']

RE_SNAKE_CASE = re.compile("[a-z0-9]+(?:_[a-z0-9]+)*")
RE_CAMEL_CASE = re.compile("([A-Z][a-z0-9]*)+")

def reportError(token, msg, errorId):
    """ Report an incorrect naming convention.
    """
    cppcheckdata.reportError(token, 'style', msg, 'naming', errorId)

def main():
    """
    Main script usable by cppcheck
    """
    parser = cppcheckdata.ArgumentParser()
    args = parser.parse_args()

    dump_files, _ = cppcheckdata.get_files(args)

    for arg in dump_files:
        if not arg.endswith('.dump') or any(f'{f}.dump' == arg for f in file_exceptions):
            continue
        print(f'Checking {arg}...')
        data = cppcheckdata.parsedump(arg)

        fname = Path(arg).stem
        res = re.match(RE_SNAKE_CASE, fname)
        if not res:
            tok = next(tok for tok in data.rawTokens)
            filename = arg.removesuffix(".dump")
            reportError(tok,
                        f'File {filename} violates snake-case naming convention', 'varname')

        # Find any use of the DDC keywords listed in ddc_keyword_map in the raw code.
        ddc_usage = [tok for tok in data.rawTokens if tok.str in ddc_keyword_map]
        for tok in ddc_usage:
            key = tok.str
            new_key = ddc_keyword_map[key]
            reportError(tok,
                        f'Please use the new naming conventions instead of DDC names ({key}->{new_key})',
                        'DDCNames')

        # Loop over configurations. A configuration is created by an ifdef
        for cfg in data.iterconfigurations():
            print(f'Checking {arg}, config {cfg.name}...')

            # Check that variables respect naming conventions
            for var in cfg.variables:
                # If anonymous variable then move on to next element
                if var.nameToken is None:
                    continue

                # Check that variables use snake case by checking against a regex
                res = re.match(RE_SNAKE_CASE, var.nameToken.str)
                if not res:
                    reportError(var.typeStartToken,
                                f'Variable {var.nameToken.str} violates snake-case naming convention',
                                'varname')

                # For class variables check that member variables start with m_
                # and that static variables start with s_
                if var.access in ('Private', 'Protected', 'Public'):
                    name = var.nameToken.str
                    if not ((name.startswith('m_') and not var.isStatic) or
                            (name.startswith('s_') and var.isStatic)):
                        if var.isStatic:
                            expected = 's_' + name
                        else:
                            expected = 'm_' + name
                        reportError(var.typeStartToken,
                                    f'Private member variable {name} violates naming convention. Expected : {expected}',
                                    'privateMemberVariable')

            # Check naming conventions for functions and classes
            for scope in cfg.scopes:
                # If scope has no body then it is irrelevant (e.g. global scope)
                if scope.bodyStart is None:
                    continue

                # Examine functions
                if scope.type == 'Function':
                    function = scope.function
                    # Naming conventions do not apply to constructors and destructors
                    # whose name cannot be chosen at will
                    if function is not None and function.type in ('Constructor', 'Destructor', 'CopyConstructor', 'MoveConstructor'):
                        continue
                    # Add an exception for gtest names
                    if scope.className == 'SetUpTestSuite' or scope.className.startswith('__'):
                        continue

                    # Check that functions use snake case by checking against a regex
                    res = re.match(RE_SNAKE_CASE, scope.className)
                    if not res:
                        reportError(scope.bodyStart,
                            f'Function {scope.className} violates snake-case naming convention',
                            'functionName')

                # Examine classes and structs
                if scope.type in ('Class', 'Struct'):
                    # Check that classes use camel case by checking against a regex
                    res = re.match(RE_CAMEL_CASE, scope.className)
                    if not res:
                        reportError(scope.bodyStart,
                            f'{scope.type} {scope.className} violates camel-case naming convention',
                            f'{scope.type.lower()}Name')

if __name__ == '__main__':
    main()
    sys.exit(cppcheckdata.EXIT_CODE)

