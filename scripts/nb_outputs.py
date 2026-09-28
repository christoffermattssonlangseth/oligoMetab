#!/usr/bin/env python3
"""Dump the text outputs (stdout, text/plain of displayed tables, error names) of an executed notebook, skipping images.

    python scripts/nb_outputs.py notebooks/analysis/analysis_spatial_metabolism.ipynb [max_chars_per_output]
"""
import json
import sys

path = sys.argv[1]
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
nb = json.load(open(path, encoding="utf-8"))
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] == "markdown":
        head = "".join(c["source"]).strip().splitlines()[0] if c["source"] else ""
        print(f"\n##### [{i}] {head}")
        continue
    outs = c.get("outputs", [])
    if not outs:
        continue
    print(f"\n----- [{i}] code cell outputs")
    for o in outs:
        if o.get("output_type") == "stream":
            print("".join(o["text"])[:limit])
        elif o.get("output_type") in ("display_data", "execute_result"):
            d = o.get("data", {})
            if "text/plain" in d and "image/png" not in d:
                t = "".join(d["text/plain"])
                if not t.startswith("<Figure") and not t.startswith("<AxesSubplot"):
                    print(t[:limit])
            elif "image/png" in d:
                print("[figure]")
        elif o.get("output_type") == "error":
            print("ERROR:", o.get("ename"), o.get("evalue"))
