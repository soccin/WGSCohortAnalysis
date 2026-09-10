# data/

`data/raw/` holds the inputs this project was run from:

- `tempoSNVManifest_<yymmdd>.csv` and `tempoSVManifest_<yymmdd>.csv`:
  Tempo manifests (`TID,NID,ProjNo,PATH,Sig`). `Sig` is the md5 of the file
  at `PATH`. Build new ones with `scan_tempo_outputs()` from the toolkit.
- optional projects workbook (`Project`, `Type`, `Ucode`, `LabSet`) used by
  `cohort.include_ucode` in `00.PARAMS.yml`.

Nothing under `data/` is written by the pipeline.

The manifests carry PHI: they name samples and the Tempo files behind them.
`.gitignore` excludes everything here except this README, and everything
under `results/` for the same reason. Do not push either to a remote.
