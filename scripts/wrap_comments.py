#!/usr/bin/env python3
import sys
import textwrap
from pathlib import Path


def wrap_comment_block(lines, prefix, width=80):
    # lines: list of strings without the prefix
    text = "\n".join(l.rstrip() for l in lines)
    # collapse multiple spaces/newlines
    text = "\n".join([t.strip() for t in text.splitlines()])
    wrapper = textwrap.TextWrapper(width=width - len(prefix), break_long_words=False, replace_whitespace=True)
    wrapped = []
    for para in text.split('\n'):
        if not para:
            wrapped.append(prefix.rstrip())
        else:
            wrapped.extend([prefix + w for w in wrapper.wrap(para)])
    return [w + "\n" for w in wrapped]


def process_r_file(path: Path):
    s = path.read_text()
    out = []
    buf = []
    buf_prefix = None
    for line in s.splitlines(keepends=True):
        if line.lstrip().startswith("#"):
            # comment or roxygen
            # capture prefix including leading whitespace and the leading #' or #
            stripped = line.lstrip()
            indent = line[:len(line) - len(stripped)]
            # find comment start (e.g., #' or #)
            if stripped.startswith("#'"):
                prefix = indent + "#' "
                content = stripped[2:].lstrip()
            else:
                prefix = indent + "# "
                content = stripped[1:].lstrip()
            if buf_prefix is None:
                buf_prefix = prefix
                buf = [content.rstrip('\n')]
            elif prefix == buf_prefix:
                buf.append(content.rstrip('\n'))
            else:
                out.extend(wrap_comment_block(buf, buf_prefix))
                buf_prefix = prefix
                buf = [content.rstrip('\n')]
        else:
            if buf:
                out.extend(wrap_comment_block(buf, buf_prefix))
                buf = []
                buf_prefix = None
            out.append(line)
    if buf:
        out.extend(wrap_comment_block(buf, buf_prefix))
    path.write_text(''.join(out))


def process_rmd_file(path: Path):
    s = path.read_text()
    out = []
    in_code = False
    buf = []
    for line in s.splitlines(keepends=True):
        if line.startswith('```'):
            if buf:
                # wrap prose buffer
                wrapped = []
                for para in '\n'.join([b.rstrip() for b in buf]).split('\n'):
                    para = para.strip()
                    if not para:
                        wrapped.append('\n')
                    else:
                        wrapped.extend([w + '\n' for w in textwrap.wrap(para, width=80, replace_whitespace=True, break_long_words=False)])
                out.extend(wrapped)
                buf = []
            out.append(line)
            in_code = not in_code
            continue
        if in_code:
            out.append(line)
        else:
            # prose
            if line.strip() == '':
                if buf:
                    buf.append('')
                else:
                    out.append(line)
            else:
                buf.append(line.rstrip('\n'))
    if buf:
        wrapped = []
        for para in '\n'.join([b.rstrip() for b in buf]).split('\n'):
            para = para.strip()
            if not para:
                wrapped.append('\n')
            else:
                wrapped.extend([w + '\n' for w in textwrap.wrap(para, width=80, replace_whitespace=True, break_long_words=False)])
        out.extend(wrapped)
    path.write_text(''.join(out))


def main(paths):
    for p in paths:
        path = Path(p)
        if not path.exists():
            continue
        if path.suffix == '.R' or path.suffix == '.Rd':
            process_r_file(path)
        elif path.suffix == '.Rmd':
            process_rmd_file(path)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        print('Usage: wrap_comments.py <file1> [file2 ...]')
