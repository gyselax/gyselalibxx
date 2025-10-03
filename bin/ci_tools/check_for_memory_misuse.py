"""
Check that CPU/GPU memory patterns are respected.
This file should be run using ./bin/run_cppcheck
"""
import sys

import cppcheckdata #pylint: disable=import-error

def reportError(token, msg, errorId, severity='error'):
    """ Report an incorrect naming convention.
    """
    cppcheckdata.reportError(token, severity, msg, 'memory', errorId)

def main():
    """
    Main script usable by cppcheck
    """
    parser = cppcheckdata.ArgumentParser()
    args = parser.parse_args()

    dump_files, _ = cppcheckdata.get_files(args)

    data_dump = {f.removesuffix('.dump'): cppcheckdata.parsedump(f) for f in dump_files}

    for f, data in data_dump.items():
        # Loop over configurations. A configuration is created by an ifdef
        for cfg in data.iterconfigurations():
            # cfg properties : clang_warnings', 'containers', 'directives', 'functions', 'macro_usage', 'name', 'preprocessor_if_conditions', 'scopes', 'setIdMap', 'set_id_map', 'set_tokens_links', 'standards', 'tokenlist', 'typedefInfo', 'valueflow', 'variables'

            class_scope = None

            # Sort scopes hiererachically
            scopes_to_add = set(cfg.scopes)
            scopes = []
            while scopes_to_add:
                valid_scopes = [s for s in scopes_to_add if s.nestedIn is None or s.nestedIn in scopes]
                scopes.extend(valid_scopes)
                scopes_to_add.difference_update(valid_scopes)

            # Find execution space for each scope
            for scope in scopes:
                # Only functions and lambdas have an execution space
                if scope.type not in ('Function', 'Unconditional'):
                    continue

                # Find body start token
                tok_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == scope.bodyStartId)
                start_idx = next((tok_idx-i for i, tok in enumerate(reversed(cfg.tokenlist[:tok_idx]),1) if tok.str in ('{', '}', ';')), 0)
                # Find any KOKKOS macros in the declaration
                exec_type = next((tok.str for tok in cfg.tokenlist[start_idx:tok_idx] if 'KOKKOS_' in tok.str), 'CPU')

                # Check parent scope. A bare for in a GPU loop is still on GPU
                parent_scope = scope
                while parent_scope.nestedIn and exec_type == 'CPU' and parent_scope.nestedIn.type in ('Function', 'Unconditional'):
                    parent_scope = parent_scope.nestedIn
                    if parent_scope.exec_space == 'GPU':
                        exec_type = parent_scope.exec_type

                # Set the exec space (CPU/GPU) and the exec_type (KOKKOS_FUNCTION, KOKKOS_LAMBDA, etc)
                scope.exec_type = exec_type
                scope.exec_space = 'CPU' if exec_type == 'CPU' else 'GPU'

                if scope.function and scope.exec_space == 'GPU' and (scope.function.hasVirtualSpecifier or scope.function.isImplicitlyVirtual):
                    tok = next(tok for tok in cfg.tokenlist if tok.Id == scope.bodyStartId).previous
                    reportError(tok,
                        f'Virtual functions cannot be reliably called from GPU',
                        'virtualFunctionOnGPU')

            # Check for bad class access from GPU
            for scope in cfg.scopes:
                if scope.type not in ('Class', 'Struct'):
                    continue

                # Find class methods
                function_scopes = [s for s in scope.nestedList if s.type == 'Function']
                functions = [s.function for s in function_scopes]

                # Check if an instance of the class can be created on GPU
                has_gpu_constructor = any(s.exec_space != 'CPU' for f,s in zip(functions, function_scopes) \
                                        if f.type == 'Constructor')

                nested_scopes = scope.nestedList
                # If any methods run on GPU and the class contains references which cannot be created
                # from GPU objects (because there is no GPU constructor)
                if any(getattr(s, 'exec_space', None) == 'GPU' for s in nested_scopes) and \
                        any(v.isReference for v in scope.varlist) and \
                        not has_gpu_constructor:
                    tok = next(tok for tok in cfg.tokenlist if tok.Id == scope.bodyStartId)
                    reportError(tok,
                        f'Class {scope.className} contains CPU references to objects.',
                        'classExecSpace')

            # Find memory space for variables
            for scope in cfg.scopes:
                for var in scope.varlist:
                    tok_start_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == var.typeStartTokenId)
                    tok_end_idx = next(i for i, tok in enumerate(cfg.tokenlist) if tok.Id == var.typeEndTokenId) + 1
                    type_descr = [tok.str for tok in cfg.tokenlist[tok_start_idx:tok_end_idx]]
                    if any(t in type_descr for t in ('FieldMem', 'VectorFieldMem', 'Field', 'VectorField', 'DiscreteToCartesian')):
                        if any(s in type_descr for s in ('host_t', 'DefaultHostExecutionSpace', 'HostSpace')):
                            var.mem_space = 'CPU'
                        else:
                            var.mem_space = 'GPU'

            # Check for memory misuse
            for scope in cfg.scopes:
                if not hasattr(scope, 'exec_space'):
                    continue

                if scope.bodyStart.file != f:
                    continue

                # Find all variables used in this scope
                used_vars = {}
                body_token = scope.bodyStart
                while body_token != scope.bodyEnd:
                    if body_token.variable:
                        used_vars[body_token.variable] = body_token
                    body_token = body_token.next

                # Find class scope
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
                        possible_funcs = [f for f in cfg.functions if f.name == func_scope.className]
                        assert len(possible_funcs) == 1
                        class_scope = scope.nestedIn
                        while class_scope.type not in ('Class', 'Struct', 'Global'):
                            class_scope = class_scope.nestedIn

                class_captured = 'CLASS' in scope.exec_type or scope.nestedIn.type in ('Class', 'Struct')

                local_access = ('Local', 'Argument')
                # Check for missing class capture
                if scope.exec_space == 'GPU' and not class_captured and class_scope and \
                        any(var.access not in local_access for var in used_vars):
                    for var, tok in used_vars.items():
                        if var.access not in local_access and var.scopeId == class_scope.Id:
                            reportError(tok,
                                f'Class variable {var.nameToken.str} used from GPU, but class not captured. Please create a proxy variable.',
                                'classVarOnGPU')

                # Check for bad class capture
                if scope.exec_space == 'GPU' and class_captured:
                    # Check for nonsensical class capture
                    if class_scope is None:
                        reportError(scope.bodyStart,
                            'Loop on GPU captures class but does not appear to be in a class.',
                            'classVarOnGPU', severity='warning')

                    # Check if *this is used in the loop
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

                    # Check for bad KOKKOS_CLASS_LAMBDA use
                    if scope.exec_type == 'KOKKOS_CLASS_LAMBDA' and \
                            not all(v in used_vars for v in class_scope.varlist) and \
                            not this_used:
                        body_token = scope.bodyStart
                        while body_token != scope.bodyEnd:
                            body_token = body_token.next
                        reportError(scope.bodyStart,
                             'Loop on GPU captures class but not all variables are used. ' +
                             'This increases CPU/GPU copies and GPU memory usage unnecessarily. ' +
                             'Please prefer proxy variables.',
                             'unnecessaryGPUclassCapture', severity='warning')

                tok = scope.bodyStart
                body = [tok]
                while tok != scope.bodyEnd:
                    tok = tok.next
                    body.append(tok)

                # Check for bad exec space use
                for var in used_vars:
                    if hasattr(var, 'mem_space'):
                        if var.mem_space != scope.exec_space:
                            # Check for variables that are actually accessed in this scope
                            usage = [tok for tok in body if tok.variable == var and tok.next.str == '(' and tok.scope == scope]
                            try:
                                usage.remove(var.typeEndToken.next)
                            except ValueError:
                                pass
                            for tok in usage:
                                msg = ("Attempted memory access from the wrong execution space. "
                                      f"The current execution space is \"{scope.exec_space}\" but the variable "
                                      f"{var.nameToken.str} is stored on a space compatible with \"{var.mem_space}\"")
                                reportError(tok, msg, 'badAccessCPUGPU')

if __name__ == '__main__':
    main()
    sys.exit(cppcheckdata.EXIT_CODE)
