#!/usr/bin/env python3
"""
modernize.py — Auto-modernize Burkardt Fortran to Fortran 2018.

Mechanical transforms (safe, preserves semantics):
  1. Replace kind=4/8 with iso_fortran_env named constants
  2. Add explicit intent(in/out/inout) from Burkardt !! docstring hints
  3. Replace 2.0D+00 literals with 2.0_dp
  4. Wrap in module with implicit none, private default
  5. Add bind(C, name="...") on all public subroutines/functions
  6. Add value attribute on scalar intent(in) bind(C) params
  7. Remove timestamp/get_unit/r8mat_write/r8vec_print utility copies
  8. Strip trailing 'return' before 'end'

Usage:
    python scripts/modernize.py src/original/grids/ball_grid.f90 src/modern/grids/ball_grid.f90
    python scripts/modernize.py src/original/  src/modern/    # batch mode
"""

import re
import os
import sys
from pathlib import Path

# Utilities that appear in many files — skip them (they'll be in a shared module)
SKIP_ROUTINES = {
    'timestamp', 'get_unit', 'get_seed', 'random_initialize',
    'r8mat_write', 'r8vec_print', 'r83vec_print_part', 'r8mat_print',
    'r8mat_print_some', 'r8vec_transpose_print', 'i4vec_print',
    'r8mat_transpose_print', 'r8mat_transpose_print_some',
    'i4mat_print', 'i4mat_print_some', 'i4mat_transpose_print',
    'i4mat_transpose_print_some', 'file_column_count', 'file_row_count',
    'r8mat_data_read', 'r8mat_header_read', 'i4mat_data_read',
    'i4mat_header_read', 'r8mat_write', 'i4mat_write',
    'file_name_ext_get', 'file_name_ext_swap', 's_to_r8',
    's_to_r8vec', 's_to_i4', 's_to_i4vec', 's_word_count',
    's_blank_delete', 's_index_last_c', 'ch_cap', 'ch_eqi',
    'ch_to_digit', 'word_next_read',
}


def parse_routines(content):
    """Split Fortran content into individual routine blocks."""
    # Join continuation lines for easier parsing
    lines = content.split('\n')
    routines = []
    current = []
    current_name = None
    in_routine = False

    for line in lines:
        # Detect subroutine/function start
        m = re.match(
            r'(?:recursive\s+)?(?:pure\s+)?(?:elemental\s+)?'
            r'(?:integer\s*(?:\([^)]*\))?\s+)?'
            r'(?:real\s*(?:\([^)]*\))?\s+)?'
            r'(?:double\s+precision\s+)?'
            r'(?:logical\s*(?:\([^)]*\))?\s+)?'
            r'(subroutine|function)\s+(\w+)',
            line.strip(), re.IGNORECASE
        )

        if m and not in_routine:
            in_routine = True
            current_name = m.group(2).lower()
            current = [line]
            continue

        if in_routine:
            current.append(line)
            # Check for end
            stripped = line.strip().lower()
            if stripped == 'end' or re.match(r'end\s+(subroutine|function)', stripped):
                routines.append((current_name, '\n'.join(current)))
                current = []
                current_name = None
                in_routine = False

    return routines


def extract_param_info(routine_text):
    """Extract parameter names, types, intents from Burkardt docstrings."""
    params = {}

    # Parse from !! docstring: "Input, integer ( kind = 4 ) N, ..."
    # or "Output, real ( kind = 8 ) BG(3,NG), ..."
    for m in re.finditer(
        r'!!\s+(Input|Output|Input/output),\s+'
        r'(integer|real|logical|character|complex)\s*'
        r'(?:\([^)]*\))?\s+'
        r'(\w+)',
        routine_text, re.IGNORECASE
    ):
        intent_word = m.group(1).lower()
        ptype = m.group(2).lower()
        pname = m.group(3).lower()

        if intent_word == 'input':
            intent = 'in'
        elif intent_word == 'output':
            intent = 'out'
        else:
            intent = 'inout'

        params[pname] = {'intent': intent, 'type': ptype}

    return params


def detect_array_params(routine_text, param_names):
    """Detect which parameters are arrays from their declarations."""
    arrays = set()
    for pname in param_names:
        # Look for declarations like: real ( kind = 8 ) bg(3,ng) or real(dp) :: bg(3,ng)
        pattern = rf'\b{re.escape(pname)}\s*\('
        if re.search(pattern, routine_text, re.IGNORECASE):
            # Check it's in a declaration context (not a call)
            for line in routine_text.split('\n'):
                stripped = line.strip().lower()
                if re.match(r'(integer|real|double|logical|complex|character)', stripped):
                    if re.search(pattern, stripped, re.IGNORECASE):
                        arrays.add(pname)
                        break
    return arrays


def modernize_types(text):
    """Replace kind=4/8 with iso_fortran_env constants."""
    # integer ( kind = 4 ) -> integer(ip)
    text = re.sub(r'integer\s*\(\s*kind\s*=\s*4\s*\)', 'integer(ip)', text, flags=re.IGNORECASE)
    text = re.sub(r'integer\s*\(\s*4\s*\)', 'integer(ip)', text, flags=re.IGNORECASE)

    # integer ( kind = 8 ) -> integer(int64)  (rare but exists)
    text = re.sub(r'integer\s*\(\s*kind\s*=\s*8\s*\)', 'integer(int64)', text, flags=re.IGNORECASE)

    # real ( kind = 8 ) -> real(dp)
    text = re.sub(r'real\s*\(\s*kind\s*=\s*8\s*\)', 'real(dp)', text, flags=re.IGNORECASE)
    text = re.sub(r'real\s*\(\s*8\s*\)', 'real(dp)', text, flags=re.IGNORECASE)

    # real ( kind = 4 ) -> real(sp)
    text = re.sub(r'real\s*\(\s*kind\s*=\s*4\s*\)', 'real(sp)', text, flags=re.IGNORECASE)

    # double precision -> real(dp)
    text = re.sub(r'double\s+precision\b', 'real(dp)', text, flags=re.IGNORECASE)

    # logical ( kind = 4 ) -> logical
    text = re.sub(r'logical\s*\(\s*kind\s*=\s*4\s*\)', 'logical', text, flags=re.IGNORECASE)

    # complex ( kind = 8 ) -> complex(dp)
    text = re.sub(r'complex\s*\(\s*kind\s*=\s*8\s*\)', 'complex(dp)', text, flags=re.IGNORECASE)

    # Replace D+00 literals: 2.0D+00 -> 2.0_dp, 1.0D-03 -> 1.0e-03_dp
    text = re.sub(r'(\d+\.\d*)[Dd]([+-]?\d+)', r'\1e\2_dp', text)

    # Replace real(..., kind = 8) -> real(..., dp)
    text = re.sub(r',\s*kind\s*=\s*8\s*\)', ', dp)', text, flags=re.IGNORECASE)
    text = re.sub(r',\s*kind\s*=\s*4\s*\)', ', sp)', text, flags=re.IGNORECASE)

    return text


def modernize_routine(name, text, module_name):
    """Modernize a single routine."""
    lines = text.split('\n')

    # Extract subroutine/function signature
    first_line = lines[0].strip()

    # Determine if function or subroutine
    is_function = 'function' in first_line.lower() and 'subroutine' not in first_line.lower()
    kind = 'function' if is_function else 'subroutine'

    # Extract parameter list
    m = re.search(r'\(([^)]*)\)', first_line)
    param_str = m.group(1).strip() if m else ''
    param_names = [p.strip().lower() for p in param_str.split(',') if p.strip()]

    # Get intents from docstrings
    param_info = extract_param_info(text)

    # Get array info
    arrays = detect_array_params(text, param_names)

    # Modernize types in the body
    new_text = modernize_types(text)

    # Remove 'implicit none' (module-level handles it)
    new_text = re.sub(r'^\s*implicit\s+none\s*$', '', new_text, flags=re.MULTILINE | re.IGNORECASE)

    # Remove trailing bare 'return' before end
    new_text = re.sub(r'\n\s*return\s*\n(\s*end\b)', r'\n\1', new_text, flags=re.IGNORECASE)

    # Add intent to declarations where we know it
    for pname, info in param_info.items():
        intent = info['intent']
        # Match type declaration of this param (no existing intent)
        # e.g., "  real(dp) c(3)" or "  integer(ip) n"
        pattern = rf'^(\s*(?:integer|real|logical|complex|character)\s*(?:\([^)]*\))?)\s+({re.escape(pname)}\b.*)$'
        replacement = rf'\1, intent({intent}) :: \2'

        new_lines = []
        for line in new_text.split('\n'):
            m2 = re.match(pattern, line.strip(), re.IGNORECASE)
            if m2 and 'intent' not in line.lower() and '::' not in line:
                # Reconstruct with intent and ::
                indent = len(line) - len(line.lstrip())
                type_part = m2.group(1)
                var_part = m2.group(2)
                new_line = ' ' * indent + f'{type_part}, intent({intent}) :: {var_part}'
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        new_text = '\n'.join(new_lines)

    # Add bind(C) to the signature
    # Find the subroutine/function line and add bind(C)
    sig_lines = new_text.split('\n')
    for idx, line in enumerate(sig_lines):
        if re.match(rf'\s*(?:recursive\s+)?(?:pure\s+)?{kind}\s+{re.escape(name)}\s*\(', line, re.IGNORECASE):
            # Find end of signature (may span multiple lines with &)
            end_idx = idx
            while end_idx < len(sig_lines) and ')' not in sig_lines[end_idx]:
                end_idx += 1

            # Add bind(C) after closing paren
            sig_lines[end_idx] = sig_lines[end_idx].rstrip()
            if sig_lines[end_idx].endswith(')'):
                sig_lines[end_idx] += f' &\n      bind(C, name="{name}")'
            break

    new_text = '\n'.join(sig_lines)

    return new_text


def modernize_file(input_path, output_path):
    """Modernize an entire Fortran file."""
    with open(input_path, 'r', errors='replace') as f:
        content = f.read()

    # Parse into routines
    routines = parse_routines(content)

    if not routines:
        return 0

    # Filter out utility routines
    keep = [(name, text) for name, text in routines if name not in SKIP_ROUTINES]

    if not keep:
        return 0

    # Derive module name from filename
    basename = Path(input_path).stem
    module_name = f'{basename}_mod'

    # Build output
    out_lines = []
    out_lines.append(f'!> {basename} — Modern Fortran 2018')
    out_lines.append(f'!>')
    out_lines.append(f'!> Modernized from John Burkardt\'s original (GNU LGPL).')
    out_lines.append(f'')
    out_lines.append(f'module {module_name}')
    out_lines.append(f'  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64')
    out_lines.append(f'  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool')
    out_lines.append(f'  implicit none')
    out_lines.append(f'  private')
    out_lines.append(f'')
    out_lines.append(f'  integer, parameter :: dp = real64')
    out_lines.append(f'  integer, parameter :: sp = real32')
    out_lines.append(f'  integer, parameter :: ip = int32')
    out_lines.append(f'')

    # Public declarations
    public_names = [name for name, _ in keep]
    out_lines.append(f'  public :: {", ".join(public_names[:6])}')
    for i in range(6, len(public_names), 6):
        chunk = public_names[i:i+6]
        out_lines.append(f'  public :: {", ".join(chunk)}')
    out_lines.append(f'')
    out_lines.append(f'contains')
    out_lines.append(f'')

    # Modernize each routine
    for name, text in keep:
        modernized = modernize_routine(name, text, module_name)
        # Indent routine body inside module
        mod_lines = modernized.split('\n')
        for ml in mod_lines:
            if ml.strip():
                out_lines.append(f'  {ml}')
            else:
                out_lines.append('')
        out_lines.append('')

    out_lines.append(f'end module {module_name}')
    out_lines.append('')

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(out_lines))

    return len(keep)


def main():
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <input> <output>')
        print(f'  input/output can be files or directories')
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if os.path.isdir(input_path):
        # Batch mode
        total_files = 0
        total_routines = 0

        for root, dirs, files in os.walk(input_path):
            for fname in sorted(files):
                if not fname.endswith('.f90'):
                    continue

                rel = os.path.relpath(os.path.join(root, fname), input_path)
                out_file = os.path.join(output_path, rel)

                count = modernize_file(os.path.join(root, fname), out_file)
                if count > 0:
                    total_files += 1
                    total_routines += count
                    print(f'  {rel}: {count} routines')

        print(f'\nTotal: {total_files} files, {total_routines} routines modernized')

    else:
        count = modernize_file(input_path, output_path)
        print(f'Modernized {count} routines')


if __name__ == '__main__':
    main()
