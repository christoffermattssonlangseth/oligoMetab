#!/usr/bin/env python3
"""Notebook hygiene: valid JSON, no hard-coded home paths in source cells, env vars used are documented in .env.example."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOME_PATH = re.compile(r"""['"](/Users/|/home/|[A-Za-z]:\\\\Users\\\\)""")
ENV_VAR = re.compile(r"""os\.(?:environ\[|getenv\()['"]([A-Z0-9_]+)['"]|resolve_path\(['"]([A-Z0-9_]+)['"]""")


def main() -> int:
    documented = {line.split("=")[0].strip() for line in (ROOT / ".env.example").read_text().splitlines() if "=" in line and not line.startswith("#")}
    documented |= {"OLIGOC4B_PUBLIC_PROCESSED_DIR"}
    problems = []
    for nb_path in sorted((ROOT / "notebooks").rglob("*.ipynb")):
        if ".ipynb_checkpoints" in nb_path.parts:
            continue
        try:
            nb = json.loads(nb_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            problems.append(f"{nb_path}: invalid JSON ({exc})"); continue
        for i, cell in enumerate(nb.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue
            src = "".join(cell.get("source", []))
            if HOME_PATH.search(src):
                problems.append(f"{nb_path}: cell {i} has a hard-coded home path")
            for m in ENV_VAR.finditer(src):
                var = m.group(1) or m.group(2)
                if var not in documented:
                    problems.append(f"{nb_path}: cell {i} uses undocumented env var {var}")
    for py in sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "notebooks" / "src").glob("*.py")):
        src = py.read_text(encoding="utf-8")
        for m in ENV_VAR.finditer(src):
            var = m.group(1) or m.group(2)
            if var not in documented:
                problems.append(f"{py}: uses undocumented env var {var}")
    for p in problems:
        print(p)
    print("OK" if not problems else f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
