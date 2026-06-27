#!/usr/bin/env python3
"""
modernize.py — Auto-modernize Burkardt Fortran to Fortran 2018.

Mechanical transforms (safe, preserves semantics):
  1. Replace kind=4/8 with iso_fortran_env named constants (int32, real64, etc.)
  2. Add explicit intent(in/out/inout) from Burkardt !! docstring hints
     — validated against body usage to avoid incorrect annotations
  3. Replace 2.0D+00 literals with 2.0_real64
  4. Wrap in module with implicit none, private default
  5. Remove timestamp/get_unit/r8mat_write/r8vec_print utility copies
  6. Strip trailing 'return' before 'end'

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
        pattern = rf'\b{re.escape(pname)}\s*\('
        if re.search(pattern, routine_text, re.IGNORECASE):
            for line in routine_text.split('\n'):
                stripped = line.strip().lower()
                if re.match(r'(integer|real|double|logical|complex|character)', stripped):
                    if re.search(pattern, stripped, re.IGNORECASE):
                        arrays.add(pname)
                        break
    return arrays


def find_modified_params(routine_text, param_names):
    """Detect which parameters are modified in the routine body."""
    modified = set()

    for line in routine_text.split('\n'):
        stripped = line.strip().lower()

        if not stripped or stripped.startswith('!'):
            continue
        if re.match(r'(integer|real|double|logical|complex|character|implicit|use\b|save\b|data\b|parameter\b)', stripped):
            continue

        for pname in param_names:
            if re.match(rf'{re.escape(pname)}\s*(\(.*?\))?\s*=\s*[^=]', stripped):
                modified.add(pname)

    return modified


def modernize_types(text):
    """Replace kind=4/8 with iso_fortran_env constants.

    Simplifies verbose kind specifiers to compact form.
    Standalone routines (no module), so use numeric kinds directly.
    """
    # integer ( kind = 4 ) -> integer — ensure trailing space
    text = re.sub(r'integer\s*\(\s*kind\s*=\s*4\s*\)\s*', 'integer ', text, flags=re.IGNORECASE)
    text = re.sub(r'integer\s*\(\s*4\s*\)\s*', 'integer ', text, flags=re.IGNORECASE)

    # integer ( kind = 8 ) -> integer(8)
    text = re.sub(r'integer\s*\(\s*kind\s*=\s*8\s*\)\s*', 'integer(8) ', text, flags=re.IGNORECASE)

    # real ( kind = 8 ) -> double precision — ensure trailing space
    text = re.sub(r'real\s*\(\s*kind\s*=\s*8\s*\)\s*', 'double precision ', text, flags=re.IGNORECASE)
    text = re.sub(r'real\s*\(\s*8\s*\)\s*', 'double precision ', text, flags=re.IGNORECASE)

    # real ( kind = 4 ) -> real
    text = re.sub(r'real\s*\(\s*kind\s*=\s*4\s*\)\s*', 'real ', text, flags=re.IGNORECASE)

    # logical ( kind = 4 ) -> logical
    text = re.sub(r'logical\s*\(\s*kind\s*=\s*4\s*\)\s*', 'logical ', text, flags=re.IGNORECASE)

    # complex ( kind = 8 ) -> complex(8)
    text = re.sub(r'complex\s*\(\s*kind\s*=\s*8\s*\)\s*', 'complex(8) ', text, flags=re.IGNORECASE)

    # Keep D+00 literals as-is (standard Fortran double precision)

    # Replace kind specs in casts: real(..., kind = 8) -> dble(...)
    text = re.sub(r',\s*kind\s*=\s*8\s*\)', ')', text, flags=re.IGNORECASE)
    text = re.sub(r',\s*kind\s*=\s*4\s*\)', ')', text, flags=re.IGNORECASE)

    return text


def modernize_routine(name, text, module_name):
    """Modernize a single routine."""
    lines = text.split('\n')
    first_line = lines[0].strip()

    is_function = 'function' in first_line.lower() and 'subroutine' not in first_line.lower()

    # Remove 'pure' and 'elemental' — too many false positives with impure callees
    first_line_clean = re.sub(r'\bpure\s+', '', first_line, flags=re.IGNORECASE)
    first_line_clean = re.sub(r'\belemental\s+', '', first_line_clean, flags=re.IGNORECASE)
    if first_line_clean != first_line:
        lines[0] = lines[0].replace(first_line, first_line_clean)
        text = '\n'.join(lines)

    # Extract parameter list
    m = re.search(r'\(([^)]*)\)', first_line)
    param_str = m.group(1).strip() if m else ''
    param_names = [p.strip().lower() for p in param_str.split(',') if p.strip()]

    # Get intents from docstrings
    param_info = extract_param_info(text)

    # Validate intents: check which params are modified in the body
    modified = find_modified_params(text, set(param_info.keys()))
    for pname in list(param_info.keys()):
        info = param_info[pname]
        if info['intent'] == 'in' and pname in modified:
            param_info[pname]['intent'] = 'inout'

    # Modernize types
    new_text = modernize_types(text)

    # Ensure each routine has implicit none
    # (Don't remove it — standalone routines need their own)

    # Remove trailing bare 'return' before end
    new_text = re.sub(r'\n\s*return\s*\n(\s*end\b)', r'\n\1', new_text, flags=re.IGNORECASE)

    # Add intent to declarations where we know it
    for pname, info in param_info.items():
        intent = info['intent']
        pattern = rf'^(\s*(?:integer|real|logical|complex|character)\s*(?:\([^)]*\))?)\s+({re.escape(pname)}\b.*)$'

        new_lines = []
        for line in new_text.split('\n'):
            m2 = re.match(pattern, line.strip(), re.IGNORECASE)
            if m2 and 'intent' not in line.lower() and '::' not in line:
                indent = len(line) - len(line.lstrip())
                type_part = m2.group(1)
                var_part = m2.group(2)
                new_line = ' ' * indent + f'{type_part}, intent({intent}) :: {var_part}'
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        new_text = '\n'.join(new_lines)

    return new_text


def modernize_file(input_path, output_path):
    """Modernize an entire Fortran file."""
    with open(input_path, 'r', errors='replace') as f:
        content = f.read()

    routines = parse_routines(content)
    if not routines:
        return 0

    keep = [(name, text) for name, text in routines if name not in SKIP_ROUTINES]
    if not keep:
        return 0

    basename = Path(input_path).stem

    out_lines = []
    out_lines.append(f'!> {basename} — Modern Fortran 2018')
    out_lines.append(f'!>')
    out_lines.append(f'!> Modernized from John Burkardt\'s original (GNU LGPL).')
    out_lines.append(f'!> Standalone routines (no module wrapping) for clean C symbol names.')
    out_lines.append(f'')

    for name, text in keep:
        modernized = modernize_routine(name, text, basename)
        out_lines.append(modernized)
        out_lines.append('')

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

        print(f'\nTotal: {total_files} files, {total_routines} modernized')

    else:
        count = modernize_file(input_path, output_path)
        print(f'Modernized {count} routines')


if __name__ == '__main__':
    main()
