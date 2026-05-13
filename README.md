# scRNA-seq Analysis Pipeline

Single-nucleus RNA sequencing (snRNA-seq) data analysis pipeline for mouse brain tissue, investigating transcriptional changes across different experimental conditions (sleep deprivation paradigm).

## Repository Structure

```
.
├── #1-QC/                              # Quality control and data preprocessing
├── #2-DR/                              # Dimensionality reduction (with integration)
├── #2.1-DR-no-integration/             # Dimensionality reduction (without integration)
├── #3-Annotation/                      # Cell type annotation (with integration)
├── #3-Annotation-no-integration/       # Cell type annotation (without integration)
├── #4-DCA_by_scDC/                     # Differential composition analysis
├── #5-DEG - no FC filtered/            # Differential expression (no fold-change filter)
├── #5-DEG-FC filtered/                 # Differential expression (fold-change filtered)
├── #6-SCENIC/                          # Gene regulatory network analysis (pySCENIC)
├── #7-subset SST/                      # SST interneuron subset analysis + DEG + GO
├── #8-SST SCENIC/                      # SST-specific SCENIC regulon analysis
├── #9-subset Inh neuron/               # Full inhibitory neuron subset analysis + DEG
├── #10-pseudobulkseq/                  # Pseudobulk sequencing & correlation analysis
├── #11-pseudotime trajectory analysis/ # Pseudotime trajectory (Slingshot & Monocle2)
├── AMP search/                         # Antimicrobial peptide gene screening
└── check TAOK1,2,3/                    # TAOK1/2/3 gene expression validation
```

## Analysis Pipeline

| Step | Module | Description |
|------|--------|-------------|
| 1 | QC | Empty droplet removal, doublet detection, mitochondrial filtering |
| 2 | Dimensionality Reduction | PCA, UMAP with and without data integration |
| 3 | Annotation | Cluster-based cell type annotation at multiple resolutions |
| 4 | DCA | Differential composition analysis across conditions |
| 5 | DEG | Differential expression with/without fold-change filtering |
| 6 | SCENIC | Transcription factor regulon inference (pySCENIC/GRNBoost2) |
| 7 | Subset SST | SST interneuron sub-clustering, DEG, GO enrichment |
| 8 | SST SCENIC | SST-specific gene regulatory network analysis |
| 9 | Subset Inh | Full inhibitory neuron subset analysis |
| 10 | Pseudobulk | Pseudobulk differential expression and correlation |
| 11 | Trajectory | Slingshot and Monocle2 pseudotime trajectory inference |

## Experimental Conditions

- **N1**: ZT6 (control)
- **N2**: SD1D (sleep deprivation 1 day)
- **N3**: SD3D (sleep deprivation 3 days)
- **N4**: SD3D R7 (sleep deprivation 3 days + recovery 7 days)
- **N5**: ZT18 (timepoint control)

## Key Tools & Packages

### R Packages
- [Seurat](https://satijalab.org/seurat/) — single-cell analysis framework
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/) — Bioconductor single-cell data container
- [scran](https://bioconductor.org/packages/scran/) — single-cell normalization and analysis
- [DropletUtils](https://bioconductor.org/packages/DropletUtils/) — droplet-based QC
- [scds](https://bioconductor.org/packages/scds/) — doublet detection
- [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) — heatmap visualization
- [Slingshot](https://bioconductor.org/packages/slingshot/) — pseudotime trajectory inference
- [Monocle2](http://cole-trapnell-lab.github.io/monocle-release/) — pseudotime analysis

### Python Tools
- [pySCENIC](https://pyscenic.readthedocs.io/) — gene regulatory network inference
- [SCENIC](https://scenic.aertslab.org/) — single-cell regulatory network inference

## Data Files

Large intermediate data files (`.rds`, `.feather`, `.loom`, `.mtx.gz`, database files) are excluded from this repository via `.gitignore`. These include:

- Seurat object files (`.rds`)
- SCENIC cisTarget databases (`.feather`)
- SCENIC expression matrices and loom files
- Raw 10x Genomics matrix files (`.mtx.gz`)

## Notes

- R scripts are stored as `.R` or `.r` files within each module folder
- Results include PDF/PNG figures and CSV/Excel summary tables
- SCENIC analyses were run on an HPC cluster with SLURM job scheduler
- The SCENIC pipeline uses the `arboreto_with_multiprocessing.py` wrapper for GRNBoost2

## References

- Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*
- Van de Sande et al. (2020) A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature Protocols*
- Street et al. (2018) Slingshot: cell lineage and pseudotime inference. *BMC Genomics*
- Qiu et al. (2017) Reversed graph embedding resolves complex single-cell trajectories. *Nature Methods*