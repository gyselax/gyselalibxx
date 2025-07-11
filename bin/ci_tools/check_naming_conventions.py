import argparse
from pathlib import Path
import re
import sys

import cppcheckdata

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

class FileErrInfo:
    def __init__(self, fname):
        self.file = fname
        self.linenr = 0
        self.column = None

def reportError(token, msg, errorId):
    cppcheckdata.reportError(token, 'style', msg, 'naming', errorId)

def main():
    parser = cppcheckdata.ArgumentParser()
    args = parser.parse_args()

    dump_files, ctu_info_files = cppcheckdata.get_files(args)

    for arg in dump_files:
        if not arg.endswith('.dump') or any(f'{f}.dump' == arg for f in file_exceptions):
            continue
        print(f'Checking {arg}...')
        data = cppcheckdata.parsedump(arg)

        # No token to report on
        #fname = Path(arg).stem
        #res = re.match(RE_SNAKE_CASE, fname)
        #if not res:
        #    filename = arg.removesuffix(".dump")
        #    reportError(token,
        #                f'File {filename} violates snake-case naming convention', 'varname')
        ddc_usage = [tok for tok in data.rawTokens if tok.str in ddc_keyword_map]
        for tok in ddc_usage:
            key = tok.str
            new_key = ddc_keyword_map[key]
            reportError(tok,
                        f'Please use the new naming conventions instead of DDC names ({key}->{new_key})',
                        'DDCNames')

        for cfg in data.iterconfigurations():
            print(f'Checking {arg}, config {cfg.name}...')


            for var in cfg.variables:
                if var.nameToken is None:
                    continue
                if var.nameToken:
                    res = re.match(RE_SNAKE_CASE, var.nameToken.str)
                    if not res:
                        reportError(var.typeStartToken,
                                    f'Variable {var.nameToken.str} violates snake-case naming convention',
                                    'varname')
                if (var.access is None) or var.access != 'Private':
                    continue

                if not ((var.nameToken.str.startswith('m_') and not var.isStatic) or
                        (var.nameToken.str.startswith('s_') and var.isStatic)):
                    reportError(var.typeStartToken,
                                f'Private member variable {var.nameToken.str} violates naming convention',
                                'privateMemberVariable')

            for scope in cfg.scopes:
                if scope.bodyStart is None:
                    continue
                if scope.type == 'Function':
                    function = scope.function
                    if function is not None and function.type in ('Constructor', 'Destructor', 'CopyConstructor', 'MoveConstructor'):
                        continue
                    if scope.className == 'SetUpTestSuite' or scope.className.startswith('__'):
                        continue
                    res = re.match(RE_SNAKE_CASE, scope.className)
                    if not res:
                        reportError(scope.bodyStart,
                            f'Function {scope.className} violates snake-case naming convention',
                            'functionName')

                if scope.type in ('Class', 'Struct'):
                    function = scope.function
                    if function is not None and function.type in ('Constructor', 'Destructor', 'CopyConstructor', 'MoveConstructor'):
                        continue
                    res = re.match(RE_CAMEL_CASE, scope.className)
                    if not res:
                        reportError(scope.bodyStart,
                            f'{scope.type} {scope.className} violates camel-case naming convention',
                            f'{scope.type.lower()}Name')

if __name__ == '__main__':
    main()
    sys.exit(cppcheckdata.EXIT_CODE)

