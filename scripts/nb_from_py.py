#!/usr/bin/env python3
"""Convert a percent-format .py (``# %%`` / ``# %% [markdown]`` cells) into an .ipynb with the given kernel."""
import re
import sys

import nbformat

src, dst, kernel = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "sc")
text = open(src, encoding="utf-8").read()
parts = re.split(r"^# %%(.*)$", text, flags=re.M)
nb = nbformat.v4.new_notebook()
nb.metadata["kernelspec"] = {"display_name": kernel, "language": "python", "name": kernel}
nb.metadata["language_info"] = {"name": "python"}
cells = []
for i in range(1, len(parts), 2):
    header, body = parts[i].strip(), parts[i + 1].strip("\n")
    if "[markdown]" in header:
        md = "\n".join(line[2:] if line.startswith("# ") else (line[1:] if line.startswith("#") else line) for line in body.splitlines())
        cells.append(nbformat.v4.new_markdown_cell(md.strip()))
    else:
        cells.append(nbformat.v4.new_code_cell(body.strip()))
nb.cells = cells
nbformat.write(nb, dst)
print(f"wrote {dst}: {len(cells)} cells")
