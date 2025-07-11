import argparse
import glob
from pathlib import Path
import sys

import cppcheckdata

def reportError(token, msg, errorId, severity='error'):
    cppcheckdata.reportError(token, 'error', msg, 'memory', errorId)
    if token.file == file:
        cppcheckdata.reportError(token, 'error', msg, 'memory', errorId)

def main():
    parser = cppcheckdata.ArgumentParser()
    args = parser.parse_args()

    dump_files, ctu_info_files = cppcheckdata.get_files(args)

    data_dump = {f.removesuffix('.dump'): cppcheckdata.parsedump(f) for f in dump_files}

    # Detect execution spaces
    for f, data in data_dump.items():
        global file
        file = Path(f)
        while len(file.suffixes) > 1:
            file = file.with_suffix('')
        file = str(file)

        for cfg in data.iterconfigurations():
            # clang_warnings', 'containers', 'directives', 'functions', 'macro_usage', 'name', 'preprocessor_if_conditions', 'scopes', 'setIdMap', 'set_id_map', 'set_tokens_links', 'standards', 'tokenlist', 'typedefInfo', 'valueflow', 'variables'

            class_scope = None

            # Find execution space for each scope
            for scope in cfg.scopes:
                if scope.type not in ('Function', 'Unconditional'):
                    continue
                if scope.nestedIn.type == 'Class':
                    class_scope = scope.nestedIn
                tok_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == scope.bodyStartId)
                start_idx = next((tok_idx-i for i, tok in enumerate(reversed(cfg.tokenlist[:tok_idx]),1) if tok.str in ('{', '}', ';')), 0)
                exec_type = next((tok.str for tok in cfg.tokenlist[start_idx:tok_idx] if 'KOKKOS_' in tok.str), 'CPU')
                scope.exec_type = exec_type
                scope.exec_space = 'CPU' if exec_type == 'CPU' else 'GPU'

            for scope in cfg.scopes:
                if scope.type not in ('Class', 'Struct'):
                    continue

                tok = next(tok for tok in cfg.tokenlist if tok.Id == scope.bodyStartId)

                function_scopes = [s for s in scope.nestedList if s.type == 'Function']
                functions = [s.function for s in function_scopes]
                has_gpu_constructor = any(s.exec_space != 'CPU' for f,s in zip(functions, function_scopes) \
                                        if f.type == 'Constructor')

                nested_scopes = scope.nestedList
                if any(getattr(s, 'exec_space', None) == 'GPU' for s in nested_scopes) and \
                        any(v.isReference for v in scope.varlist) and \
                        not has_gpu_constructor:
                    reportError(tok,
                        f'Class {scope.className} contains CPU references to objects.',
                        'classExecSpace')

                for var in scope.varlist:
                    tok_start_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == var.typeStartTokenId)
                    tok_end_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == var.typeEndTokenId) + 1
                    type_descr = [tok.str for tok in cfg.tokenlist[tok_start_idx:tok_end_idx]]
                    if any(t in type_descr for t in ('FieldMem', 'VectorFieldMem', 'Field', 'VectorField', 'DiscreteToCartesian')):
                        if (s in type_descr for s in ('host_t', 'DefaultHostExecutionSpace', 'HostSpace')):
                            var.mem_space = 'CPU'
                        elif (s in type_descr for s in ('DefaultExecutionSpace',)):
                            var.mem_space = 'GPU'
                        else:
                            # Look for mem_space
                            print(type_descr)
                            pass

            for scope in cfg.scopes:
                if not hasattr(scope, 'exec_space'):
                    continue

                used_vars = {}
                body_token = scope.bodyStart
                while body_token != scope.bodyEnd:
                    if body_token.variable:
                        used_vars[body_token.variable] = body_token
                    body_token = body_token.next

                class_scope = scope.nestedIn
                while class_scope.type not in ('Class', 'Struct', 'Global'):
                    class_scope = class_scope.nestedIn

                if class_scope.type == 'Global':
                    func_scope = scope.nestedIn
                    while func_scope.type not in ('Global', 'Function'):
                        func_scope = func_scope.nestedIn
                    if func_scope.type == 'Global':
                        class_scope = None
                    else:
                        func_id = func_scope.bodyStart.previous.link.previous.functionId
                        func = next(f for f in cfg.functions if f.Id == func_id)
                        class_scope = func.nestedIn

                class_captured = 'CLASS' in scope.exec_type or scope.nestedIn.type in ('Class', 'Struct')

                local_access = ('Local', 'Argument')
                if scope.exec_space == 'GPU' and not class_captured and class_scope and \
                        any(var.access not in local_access for var in used_vars):
                    for var, tok in used_vars.items():
                        if var.access not in local_access and var.scopeId == class_scope.Id:
                            reportError(tok,
                                f'Class variable {var.nameToken.str} used from GPU, but class not captured. Please create a proxy variable.',
                                'classVarOnGPU')

                if scope.exec_space == 'GPU' and class_captured:
                    if class_scope is None:
                        reportError(scope.bodyStart,
                            f'Loop on GPU captures class but does not appear to be in a class.',
                            'classVarOnGPU', severity='warning')

                    this_used = False
                    body_token = scope.bodyStart
                    while body_token != scope.bodyEnd:
                        if body_token.valueType:
                            class_scope_id = body_token.valueType.typeScopeId
                            if class_scope_id:
                                used_class_scope = next(s for s in cfg.scopes if s.Id == class_scope_id)
                                if used_class_scope == class_scope:
                                    this_used = True
                                    break
                        if body_token.functionId:
                            func = next(f for f in cfg.functions if f.Id == body_token.functionId)
                            if func.nestedIn == class_scope:
                                this_used = True
                                break
                        body_token = body_token.next

                    if scope.exec_type == 'KOKKOS_CLASS_LAMBDA' and \
                            not all(v in used_vars for v in class_scope.varlist) and \
                            not this_used:
                        body_token = scope.bodyStart
                        while body_token != scope.bodyEnd:
                            body_token = body_token.next
                        reportError(scope.bodyStart,
                            f'Loop on GPU captures class but not all variables are used. ' +
                             'This increases CPU/GPU copies and GPU memory usage unnecessarily. ' +
                             'Please prefer proxy variables.',
                             'unnecessaryGPUclassCapture', severity='warning')

                for var in used_vars:
                    if hasattr(var, 'mem_space'):
                        if var.mem_space != scope.exec_space:
                            msg = ("Attempted memory access from the wrong execution space. "
                                  f"The current execution space is \"{scope.exec_space}\" but the variable "
                                  f"{var.nameToken.str} is stored on a space compatible with \"{v.mem_space}\"")
                            reportError(scope.bodyStart, msg, 'badAccessCPUGPU')

if __name__ == '__main__':
    main()
    sys.exit(cppcheckdata.EXIT_CODE)
