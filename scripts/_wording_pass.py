"""One-off wording pass (2026-09-30): replace group-specific phrasing with project-neutral wording in docs, sources and notebook markdown. Kept for re-application after re-executing notebooks from older sources."""
import glob, json, re, sys
REPL = [
    ("The collaborators' central question is ketone-body metabolism in the EAE spinal cord.", "The central question of the project is ketone-body metabolism in the EAE spinal cord."),
    ("**What this means for the collaborators' question.**", "**What this means for the central question.**"),
    ("| EAE ketone notebook: the collaborators' question, asked of Xenium EAE", "| EAE ketone notebook: the central question, asked of Xenium EAE"),
    ("The collaborators' central question, asked of the in-situ data alone.", "The central question of the project, asked of the in-situ data alone."),
    ("The collaborators' central question, answered with", "The central question of the project, answered with"),
    ("| The collaborators' question: ketone handling in Xenium EAE", "| The central question: ketone handling in Xenium EAE"),
    ("## 6c. The collaborators' question: ketone handling", "## 6c. The central question: ketone handling"),
    ("**Implication for the collaborators.**", "**Implication.**"),
    ("Five answers for the collaborators", "Five answers"),
    ("collaborators'", "project's"), ("collaborators", "the project"),
]
def fix(t):
    for a, b in REPL:
        t = t.replace(a, b)
    return t
n = 0
for f in glob.glob("docs/*.md") + ["README.md", "report/report.md"] + glob.glob("notebooks/src/*.py"):
    s = open(f).read(); t = fix(s)
    if t != s: open(f, "w").write(t); n += 1
for f in glob.glob("notebooks/analysis/*.ipynb"):
    nb = json.load(open(f)); ch = False
    for c in nb["cells"]:
        if c["cell_type"] == "markdown":
            src = "".join(c["source"]); t = fix(src)
            if t != src: c["source"] = t; ch = True
    if ch: json.dump(nb, open(f, "w"), indent=1, ensure_ascii=False); n += 1
print("files changed:", n)
import subprocess; print("remaining:", subprocess.run("grep -rli collaborator --include='*.md' --include='*.py' --include='*.ipynb' . | grep -v '^./.git' | wc -l", shell=True, capture_output=True, text=True).stdout.strip())
