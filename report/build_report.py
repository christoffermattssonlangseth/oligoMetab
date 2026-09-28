#!/usr/bin/env python3
"""Render report/report.md to report/OligoMetab_metabolism_report.pdf via HTML and headless Chrome.

    python report/collect_figures.py      # pull figures from the executed notebooks
    python report/build_report.py         # write report.html and the PDF

Requires the `markdown` Python package and Google Chrome (macOS path below; override with CHROME=/path/to/chrome).
"""
from __future__ import annotations

import os
import subprocess
from datetime import date
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent
MD = ROOT / "report.md"
HTML = ROOT / "report.html"
PDF = ROOT / "OligoMetab_metabolism_report.pdf"
CHROME = os.getenv("CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")

CSS = """
@page { size: A4; margin: 18mm 16mm 18mm 16mm; }
body { font-family: -apple-system, "Helvetica Neue", Helvetica, Arial, sans-serif; font-size: 10.5pt; line-height: 1.45; color: #1a1a1a; max-width: 178mm; margin: 0 auto; }
h1 { font-size: 20pt; margin: 0 0 4pt 0; }
h2 { font-size: 14pt; margin: 22pt 0 6pt 0; border-bottom: 1px solid #bbb; padding-bottom: 2pt; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 14pt 0 4pt 0; page-break-after: avoid; }
p { margin: 5pt 0; }
ul, ol { margin: 4pt 0 6pt 0; padding-left: 18pt; }
li { margin: 2pt 0; }
table { border-collapse: collapse; font-size: 9pt; margin: 6pt 0 10pt 0; width: 100%; page-break-inside: avoid; }
th, td { border: 1px solid #ccc; padding: 3pt 5pt; text-align: left; vertical-align: top; }
th { background: #f0f0f0; }
figure { margin: 10pt 0 12pt 0; page-break-inside: avoid; text-align: center; }
figure img { max-width: 100%; max-height: 118mm; }
figcaption { font-size: 9pt; color: #444; text-align: left; margin-top: 4pt; }
.subtitle { color: #555; font-size: 11pt; margin-bottom: 14pt; }
.box { background: #f6f7f9; border-left: 3px solid #4a6fa5; padding: 8pt 10pt; margin: 8pt 0 12pt 0; }
.pagebreak { page-break-before: always; }
code { font-size: 9pt; background: #f3f3f3; padding: 0 2pt; }
"""


def main() -> int:
    text = MD.read_text(encoding="utf-8").replace("{{DATE}}", date.today().isoformat())
    body = markdown.markdown(text, extensions=["tables", "attr_list", "md_in_html", "fenced_code"])
    html = f"<!doctype html><html><head><meta charset='utf-8'><title>OligoMetab metabolism report</title><style>{CSS}</style></head><body>{body}</body></html>"
    HTML.write_text(html, encoding="utf-8")
    cmd = [CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={PDF}", HTML.resolve().as_uri()]
    subprocess.run(cmd, check=True, capture_output=True, timeout=180)
    print(f"wrote {PDF} ({PDF.stat().st_size/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
