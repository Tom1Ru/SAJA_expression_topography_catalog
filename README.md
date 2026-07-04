# Expression topography catalog

A4 PDF catalog (0.5 page per gene) with 8×8 expression topography, genotype 5×5 panels, LM parameters, and biomarker metrics.

[expression_topography_catalog_PR1_VSP1.pdf](output/figures/expression_topography_catalog_PR1_VSP1.pdf)

## Layout (per gene)

- Header band: Transcript ID, GeneSymbol, ShortDescription, CuratorSummary
- Top-left: 8×8 mean expression contour (log2 RPM+1, 64 conditions)
- Bottom-left: 5×5 genotype heatmaps (Col-0, npr3-2 npr4-2, npr3-1 npr4-3)

- Cluster: Assigned gene cluster name
- SA_r: Correlation coefficient with SA concentration
- JA_r: Correlation coefficient with JA concentration

- SA_Sel_freq: Selection frequency as an SA biomarker
- JA_Sel_freq: Selection frequency as a JA biomarker

- SA_peak: SA concentration at which expression reaches its peak
- JA_peak: JA concentration at which expression reaches its peak

- Category: Classification of the SA/JA response topography
- r2: Coefficient of determination of the linear model
- adj_r2: Adjusted coefficient of determination of the linear model

- coef_intercept: Intercept of the linear model
- coef_SA: Coefficient for the main effect of SA
- coef_JA: Coefficient for the main effect of JA
- coef_SA:JA: Coefficient for the SA × JA interaction term

## Requirements

R packages: `ggplot2`, `grid`, `patchwork`, `viridis`, `scales`, `dplyr` .

PDFs are written to `output/figures/`.

## PR1 + VSP1 catalog

```bash
cd expression_topography_catalog
Rscript run_PR1_VSP1_catalog.R
```

Output: [output/figures/expression_topography_catalog_PR1_VSP1.pdf](output/figures/expression_topography_catalog_PR1_VSP1.pdf)

## Test (PR1, VSP1, JAZ1, JAZ2, JAZ6)

```bash
cd expression_topography_catalog
Rscript run_test_catalog.R
```

Output: `output/figures/expression_topography_catalog_test_PR1_VSP1_JAZ.pdf` 

## Full SA/JA responsive genes catalog (4,939 genes)


```bash
Rscript run_DEGn3_full_catalog.R
```

Output: `output/figures/expression_topography_catalog_DEGn3_4939genes.pdf`

## Custom gene list

```bash
Rscript run_catalog.R my_genes.txt output/figures/my_catalog.pdf
```

`my_genes.txt`: one AGI ID per line (e.g. `AT2G14610`).

## Data (`data/`)

| File | Content |
|------|---------|
| `gene_annotation.csv` | Araport11 annotation |
| `gene_lm_stats.csv` | 8×8 LM classification |
| `gene_cluster_assignments.csv` | DEGn3 55 micro clusters |
| `gene_metrics.csv` | SA/JA r, select freq, SA/JA peak (mM) |
| `mean64_log2rpm.rds` | Expressed genes × 64 conditions |
| `genotype_5x5_means.csv` | 5×5 means per genotype (6 hr) |

## Repository layout

```
config.R                  # paths and display constants
run_PR1_VSP1_catalog.R    # PR1 + VSP1 PDF (linked above)
run_test_catalog.R        # small test PDF
run_DEGn3_full_catalog.R  # full DEGn3 PDF
run_catalog.R             # custom gene list → PDF
R/
  make_catalog_pdf.R      # PDF assembly
  plot_catalog_panel.R    # single-gene panel
  contour_utils.R         # 8×8 contour plots
  data_loaders.R          # load data/
  vec64_utils.R           # 64-condition vectors
  order_genes_by_cluster.R
data/                     # pre-built inputs (required)
output/figures/           # generated PDFs (gitignored except PR1_VSP1 sample)
```
