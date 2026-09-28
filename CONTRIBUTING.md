# Contributing

Same conventions as `OligoC4b`:

- Notebook-first, but repeated logic lives in `scripts/oligometab.py` (gene panel, helpers). Notebooks read data from the env vars in `.env.example` or from `../../data/`.
- Notebook sources are kept as percent-format Python in `notebooks/src/` and converted with `python scripts/nb_from_py.py <src.py> <out.ipynb> sc`; edit the source, regenerate, then execute with `scripts/run_notebooks.sh` (kernel `sc`). Commit executed notebooks *with* outputs.
- Summary tables written by the notebooks go to `results/` (ignored by git except for the small CSVs listed in `.gitignore`).
- No raw data, secrets or machine-specific absolute paths in source cells; run `python scripts/check_notebooks.py` before committing (also runs in CI).
- When adding a notebook: row in `README.md`, summary in `docs/notebooks.md`; findings go to `docs/metabolism_findings.md`.
