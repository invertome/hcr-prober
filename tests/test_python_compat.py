"""Static-analysis regressions for Python version compatibility.

PEP 701 (Python 3.12+) lifted the restriction that f-string expressions
cannot reuse the same quote character as the surrounding string. Code
that does this (e.g. f'... {", ".join(x)} ...' with both the f-string
and the inner string using ' or both using ") parses fine on 3.12+
but fails with SyntaxError on 3.8-3.11.

setup.py declares python_requires>=3.8, so we must keep the source
parseable on 3.8-3.11. This test scans every shipped .py file for the
offending pattern and fails if any are reintroduced.
"""
import pathlib
import re

PKG = pathlib.Path(__file__).resolve().parent.parent / 'hcr_prober'

# Same-quote-nested f-string. Greedy on inner content, but rejects '}'
# in the expression body so we don't span across format-spec boundaries.
SAME_QUOTE_FSTRING = [
    re.compile(r"f'[^']*\{[^}]*'[^}]*'[^}]*\}"),
    re.compile(r'f"[^"]*\{[^}]*"[^}]*"[^}]*\}'),
]


def test_no_same_quote_nested_fstrings():
    offenders = []
    for path in PKG.rglob('*.py'):
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            for pat in SAME_QUOTE_FSTRING:
                if pat.search(line):
                    offenders.append((str(path.relative_to(PKG.parent)), lineno, line.strip()))
                    break
    if offenders:
        msg_lines = [
            'Same-quote-nested f-strings (Python 3.12+ only) found:',
            *[f'  {p}:{n}  {ln}' for p, n, ln in offenders],
            '',
            'Fix by changing the outer f-string quote, escaping, or pre-computing.',
            'Examples:',
            '  BAD : f\'list: {\\\', \\\'.join(items)}\'',
            '  GOOD: f"list: {\', \'.join(items)}"',
        ]
        raise AssertionError('\n'.join(msg_lines))
