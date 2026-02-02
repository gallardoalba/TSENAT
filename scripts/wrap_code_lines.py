#!/usr/bin/env python3
"""Auto-wrap long code lines in R tests and source files.

This script targets files under `tests/`, `R/`, and `vignettes/` and attempts a
conservative wrap for lines longer than `width` by splitting argument lists at
commas. It's designed for test files and should be reviewed before committing.
"""
import io
import os
import re
import sys

width = 80
paths = ["tests", "R", "vignettes"]

def wrap_line(line, indent="  "):
    if len(line) <= width:
        return line
    if "(" not in line or "," not in line:
        return line
    # naive: split at commas, but avoid splitting inside quotes
    parts = []
    cur = []
    in_sq = False
    in_dq = False
    esc = False
    for ch in line:
        if esc:
            cur.append(ch)
            esc = False
            continue
        if ch == "\\":
            cur.append(ch)
            esc = True
            continue
        if ch == "'" and not in_dq:
            in_sq = not in_sq
            cur.append(ch)
            continue
        if ch == '"' and not in_sq:
            in_dq = not in_dq
            cur.append(ch)
            continue
        if ch == "," and not in_sq and not in_dq:
            cur.append(ch)
            parts.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append(''.join(cur))

    # If only one part, nothing to do
    if len(parts) <= 1:
        return line

    # Find first non-space prefix (indent)
    m = re.match(r"^(\s*)", line)
    base_indent = m.group(1) if m else ""

    # Rebuild: put first part on same line, remaining on new lines with extra indent
    first = parts[0].rstrip('\n')
    rest = [p.strip() for p in parts[1:]]
    wrapped = first + '\n'
    for i, p in enumerate(rest):
        sep = ''
        # keep trailing characters (closing parens, etc.) on last part
        wrapped += base_indent + indent + p
        if i < len(rest) - 1:
            wrapped += '\n'
    # Ensure final newline
    if not wrapped.endswith('\n'):
        wrapped += '\n'
    return wrapped


def process_file(path):
    try:
        with io.open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception:
        return False, []

    changed = False
    out = []
    for line in lines:
        if len(line.rstrip('\n')) > width:
            new = wrap_line(line)
            if new != line:
                changed = True
                out.append(new)
            else:
                out.append(line)
        else:
            out.append(line)

    if changed:
        # backup
        try:
            bkp = path + '.bak'
            if not os.path.exists(bkp):
                with io.open(bkp, 'w', encoding='utf-8') as fb:
                    fb.writelines(lines)
        except Exception:
            pass
        with io.open(path, 'w', encoding='utf-8') as f:
            f.writelines(out)
    return changed, out


def main():
    changed_files = []
    for base in paths:
        if not os.path.isdir(base):
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            # skip Rcheck directories
            if '.Rcheck' in dirpath or '..Rcheck' in dirpath:
                continue
            for fn in filenames:
                if not fn.endswith(('.R', '.r', '.Rmd')):
                    continue
                path = os.path.join(dirpath, fn)
                ok, _ = process_file(path)
                if ok:
                    changed_files.append(path)
    if changed_files:
        print('Wrapped lines in:')
        for p in changed_files:
            print(p)
        return 0
    print('No changes')
    return 0


if __name__ == '__main__':
    sys.exit(main())
