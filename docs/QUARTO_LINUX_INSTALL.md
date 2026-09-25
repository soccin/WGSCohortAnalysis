# Quarto on RHEL 8 without root

For a future HTML report (`templates/project/analysis/report.qmd`).
Nothing in the toolkit needs quarto today; the pipeline writes xlsx and
PDF.

## 1. Install the CLI tarball in user space

```
mkdir -p ~/opt ~/bin
cd ~/opt
# pick the current release from https://github.com/quarto-dev/quarto-cli/releases
VER=1.6.42
curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v${VER}/quarto-${VER}-linux-amd64.tar.gz
tar xzf quarto-${VER}-linux-amd64.tar.gz
ln -sfn ~/opt/quarto-${VER} ~/opt/quarto
ln -sf ~/opt/quarto/bin/quarto ~/bin/quarto
```

The tarball is self-contained (bundles Deno and its own pandoc). On a
firewalled node download it elsewhere and copy it over.

## 2. PATH and modules

Add to `~/.bashrc` (or an environment module):

```
export PATH=$HOME/bin:$PATH
```

If R is loaded through `module load R/4.5.1`, do that before rendering so
quarto finds the same `Rscript`. Quarto uses `QUARTO_R` when set:

```
export QUARTO_R=$(dirname $(which Rscript))
```

## 3. Verify

```
quarto check
```

It reports the pandoc it bundles; the existing `~/bin/pandoc` is not used
unless `QUARTO_PANDOC` points at it. Either is fine. `quarto check` also
confirms R, `rmarkdown` and `knitr` are visible.

## 4. R packages

```r
install.packages(c("quarto", "rmarkdown", "knitr", "DT", "gt"))
```

`quarto::quarto_render("report.qmd")` renders from R; the CLI does the same.

## 5. Shape of the future report.qmd

A project-level `analysis/report.qmd` would read `00.PARAMS.yml`, load the
toolkit, and pull the stage outputs rather than recompute. Report sources
(`.qmd`, `.Rmd`) go in `analysis/`, not the project root. knitr runs the
code with `analysis/` as the working directory, so paths are built from
`here::here()`, which finds the root through the template's `.here` file:

````
---
title: "`r PARAMS$prefix` cohort report"
format:
  html:
    embed-resources: true
    toc: true
---

```{r setup}
root <- here::here()
source(file.path(Sys.getenv("WCA_HOME"), "load.R"))
PARAMS <- wca_read_params(fs::path(root, "00.PARAMS.yml"))
s <- function(x) stage_load(PARAMS, "04_summaries", x)
```

## Cohort
```{r}
s("samples") |> DT::datatable()
```

## Recurrent genes
```{r}
plot_recurrent_genes(s("summary_all"))
```

## Oncoprint
```{r fig.height=10}
plot_oncoprint(stage_load(PARAMS, "03_events", "events"))
```
````

`embed-resources: true` gives a single self-contained HTML file that can be
copied off the cluster. Render from inside `analysis/`, into the run's
results folder:

```
cd analysis
quarto render report.qmd --output-dir ../results/run01 \
  --output <prefix>_Report_<yymmdd>.html
```

Rendering `analysis/report.qmd` from the project root with `--output-dir`
fails in Quarto 1.10.18: the chunks run, then post-processing looks for
`report_files/libs/quarto-html/quarto.js` in the directory it was started
from instead of beside the qmd, and leaves an incomplete html in the root.

Add `scripts/07_report_html.R` calling `quarto::quarto_render()` when this
is wanted. Untested: whether `quarto_render()` hits the same failure when
called from the project root; it runs the same CLI.
