# AGENTS.md

Julia ModelingToolkit (MTK) model of CaMKII signaling in neonatal rat ventricular myocytes
(NRVM), plus Literate notebooks published as a Jupyter Book site (GitHub Pages).
This is a scientific artifact repo: there is **no test suite**; executing the `docs/`
notebooks in CI is the verification.

## Toolchain

- Julia 1.13 (`compat` in both Project.tomls); CI pins the exact patch from `julia_version`
  in `Manifest.toml` (currently 1.13.2).
- Root is a Pkg **workspace** (`[workspace] projects = ["Model"]`, `[sources] Model = {path="Model"}`).
  Always run Julia with `--project=@.` from the repo root; the `Model` package in `Model/`
  resolves through the workspace.
- Setup: `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`
- To manually update Julia packages: `julia --project=. -e 'using Pkg; Pkg.update()'`
- `Manifest.toml` is committed; a weekly bot opens "Julia Dependency Update" PRs (automerge
  label, squash-merged via Kodiak). Don't hand-edit dependency versions.

## Notebooks: the `docs/**/*.jl` files are the source of truth

- `docs/**/*.jl` are Literate scripts. `docs/**/*.ipynb` are generated, untracked build
  artifacts — never edit them. Worse, `.github/ci.jl` regenerates the `.jl` **from** any stale
  `.ipynb` whose SHA misses the cache, so a leftover `.ipynb` in `docs/` overwrites the source script.
- Run everything the way CI does (parallel, results cached by SHA256 in `.cache/`):
  `julia --project=@. -p 8 .github/ci.jl` — from repo root; env `NBCACHE` (default `.cache`),
  `DOCDIR` (default `docs`). CI then syncs with `cp -rf .cache/docs/* docs/`.
  `.cache/` and generated `docs/*.ipynb` are gitignored artifacts.
- A notebook is skipped when `.cache/docs/<name>.sha` matches the script's SHA256 — delete
  that `.sha` (and cached `.ipynb`) to force a re-run after behavior-preserving edits.
- Run one notebook without the cache machinery:
  `julia --project=@. -e 'ENV["GKSwstype"]="100"; using Literate; Literate.notebook("docs/fig04-freq.jl", "docs"; mdstrings=true, execute=true)'`
  (`GKSwstype=100` for headless GR plots; the output-dir arg must match the file's dir so the
  `.ipynb` lands next to the `.jl`).
- Site build (needs Node): `cd docs && jupyter-book build --html` → `docs/_build/html`.
  `docs/myst.yml` TOC pulls in `*.ipynb` and `fits/*.ipynb`, so notebooks must be executed first.
  CI sets `BASE_URL=/<repo-name>` for Pages.
- Heavy notebooks: `using Model` builds and `mtkcompile`s the whole ECC system at load
  (`Model.DEFAULT_SYS`, also warmed in `Model/src/precompile.jl`); notebooks reuse it rather
  than rebuilding — do the same.

## Model package conventions (Model/src)

- Time base is milliseconds: `const second = 1000`; `μM`/`mM`/`mV` etc. are scaled constants
  (`Model/src/utils.jl`). Always write literals with unit constants (`100.0second`,
  `1.796mM`, `0.1μM`); bare SI numbers are wrong.
- Only `build_neonatal_ecc_sys` and `build_stim_callbacks` are exported. Notebooks import
  internals explicitly, e.g. `using Model: second, Hz, μM, hil`, and call builders qualified,
  e.g. `Model.get_camkii_sys(...)`.
- Current/subsystem builders follow `get_*_eqs` / `get_*_sys` returning MTK equations/`System`s,
  composed into the full model in `ecc_neonatal.jl`.
- Stiff ODEs: solve with `KenCarp47()` or `FBDF()`, tolerances 1e-8–1e-10 (match neighboring
  notebooks); steady states via `SteadyStateDiffEq`'s `DynamicSS` — see `docs/u0.jl`, which
  prints initial-condition literals.
- Fit notebooks in `docs/fits/` produce parameters baked into `Model/src`; refitting there
  implies updating model constants.

## extras/ and data (not run by CI, not in site TOC)

- `extras/ecc/*.jl`: standalone MTK snippets for currents NOT wired into the NRVM model
  (I_tof, I_Kur, I_Kr, I_Kp, I_K1, I_Cl, I_PCA, Markov I_Na/I_Ca). Reference material only.
- `extras/normalize-data.jl` regenerates `docs/data/*-normalized.csv` from the raw CSVs and
  must run from the repo root (it uses `pwd()`, not `@__DIR__`). Never hand-edit the
  `-normalized.csv` files.
- Model provenance (Chang 2019 CaMKII, Christensen ROS, Moritti ISO, Korhonen NRVM) is
  documented in `docs/index.md` and `docs/modeldesc.md`; keep the two in sync when the model
  changes. A scheduled lychee job files broken links found in `docs/**/*.md` and `docs/**/*.jl`
  as an issue, so keep links valid in both.
