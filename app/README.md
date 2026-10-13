Historical 2026-10-11 snapshot. For current manuscript coverage, see the repository root README and docs/REPRODUCING_CURRENT.md.

# Static research explorer source

This is the existing Sesamin Atlas v1.1.3 frontend and previously published derived display snapshot. It supports local viewing with:

```sh
python -m http.server 8000 --directory app
```

Open the local server in a browser. There is no build step, package installation, backend, login, new statistical analysis or automatic remote publishing. HTML/CSS/JavaScript and 14 derived data files are included. Original source data and scientific limitations are described in the app and top-level README.

The snapshot predates the manuscript's matched-scale AREG update. Its historical pooled-UMI summaries must not be used to override the newer same-scale result in `analysis/matched_scale/`. The original public developer attribution is retained. Company image assets and contact email addresses were omitted for this source release; the favicon uses the existing generic SVG. No other scientific display values were changed.

Data builders and browser/capture source are archived under `historical_pipeline/app_tools/`; they reference the earlier project layout and are not newly executed here. Remote identity/publishing scripts, access tokens, deployment logs, private visual reviews and old UI-edit scripts are excluded. No new source licence has been assigned.
