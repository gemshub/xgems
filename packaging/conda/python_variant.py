"""Print a conda-build --variants string selecting one python out of the conda-forge pinning.

The pinning zips `python` with `is_python_min`, so both must be picked at the same index.
Usage: python_variant.py <conda_build_config.yaml of conda-forge-pinning> <python version, e.g. 3.13>
Exit status 1 if that python is not (yet) in the pinning.
"""
import re
import sys


def block(text, key):
    """Return the list items (selectors stripped) of the top-level `key:` list."""
    m = re.search(rf"^{key}:\s*\n((?:[ \t]+.*\n|\n|#.*\n)*)", text, re.M)
    items = []
    for line in (m.group(1) if m else "").splitlines():
        line = line.split("#")[0].rstrip()
        if line.lstrip().startswith("- "):
            items.append(line.lstrip()[2:].strip().strip("'\""))
    return items


pinning, version = sys.argv[1], sys.argv[2]
text = open(pinning).read()
pythons, is_min = block(text, "python"), block(text, "is_python_min")
for i, spec in enumerate(pythons):
    if spec.startswith(version + ".") and i < len(is_min):
        print("{python: ['%s'], is_python_min: ['%s']}" % (spec, is_min[i]))
        sys.exit(0)
sys.exit(1)
