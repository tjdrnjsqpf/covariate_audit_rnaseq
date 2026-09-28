# covariate_audit_rnaseq

Code and summary tables for

> Lee S-G, Park C. **Sex adjustment in 2,000 public RNA-seq studies: free when balanced, costly when confounded, and rarely performed.** *Genome Biology* (submitted, 2026).

The study reanalyses 385 human and 1,483 mouse two-group bulk RNA-seq studies from recount3 with and without a sex term, compares each sex-adjusted fit with two negative controls (the sex label permuted across all samples, and permuted within each group so that the design's collinearity is preserved), induces imbalance experimentally within 85 balanced studies, repeats the permutations to ask which genes adjustment removes, and re-reads the source publications to ask whether the original analyses adjusted.

## Layout

| Path | Contents |
|---|---|
| `scripts/00_pilot/` | recount3 metadata download and parsing; sex inference from expression; the per-study differential-expression driver (`03_pilot_de_one_project.R`, `03m_pilot_de_mouse.R`: models M0–M7 and their permuted controls) |
| `scripts/01_survey/` | design packets and the language-model design classifier with its validation |
| `scripts/02_reanalysis/` | the main reanalysis runners, summary collectors, pathway analysis (fgsea), semi-synthetic simulation, release configuration |
| `scripts/03_qc/` | pseudoreplication audit, chance-imbalance null, sex-discordance sensitivity, batch survey, Δ regression |
| `scripts/04_gate2/` | publication retrieval and the language-model extraction of adjustment practice |
| `scripts/05_report/` | figure scripts (`26_draft_figures_v3.R`, `31_supp_figures_v3.R`, `32_figures_v4.R`), source-data export, manuscript builder |
| `scripts/06_revision/` | within-group control (43, 45), induced-imbalance experiment (40–42, 44, 50), repeated permutations and content of the change (46–49, 51), tissue heterogeneity (48), GTEx enrichment (52–53), and the smaller revision analyses (30–39) |
| `config/` | the release configuration defining every comparison as run (`release_{human,mouse}_runs.tsv`), the analysis sets, synthetic-experiment study lists, sex-linked gene panel, exclusions, and the adjudicated design labels |
| `results/summary/` | one row per study: sample composition, imbalance class, and every change metric for every model pair (`main_{human,mouse}_summary_v2.tsv`); project-level metadata summaries |
| `results/revision/` | class- and bin-level tables behind the Results (within-group control `W_*`, induced imbalance `S2_*`, content of the change `C_*` with per-permutation files in `content/`, GTEx `G_*`, tissue heterogeneity `T_*`, revision analyses `E*`) |
| `results/c2/` | pathway-level summaries (Hallmark and Reactome) per study |
| `results/qc/` | pseudoreplication, chance-imbalance null, discordance sensitivity, batch imbalance |
| `results/design_llm/` | design-classifier prompt, all verdicts with rationales, validation labels and review sheets |
| `results/gate2/` | publication-audit prompt, per-record extractions with evidence sentences, merged calls, verification sheet |
| `results/figures/` | figure panels (PDF/PNG) and per-panel source data for the main and supplementary figures |
| `results/simulation/` | semi-synthetic simulation outputs |
| `data/annot/` | GENCODE v29 gene-to-chromosome table used for the sex-chromosome analyses |

Not included here: the recount3 count cache (downloaded by the scripts), the per-study differential-expression tables (about 6 GB; deposited on Zenodo, see below), and the GTEx v8 sex-biased gene statistics (Oliva et al. 2020; download `GTEx_Analysis_v8_sbgenes.tar.gz` from the GTEx portal into `data/external/gtex_v8_sbgenes/`).

## Requirements

- R 4.5 with `recount3` (1.20), `edgeR` (4.8), `limma` (3.66), `sva` (3.58), `fgsea` (1.36), `msigdbr` (26.1), `data.table`, `ggplot2`, `patchwork`
- Python 3.12 with `anthropic`, `pydantic`, `requests` — only for the two language-model steps; every verdict those steps produced is included under `results/design_llm/` and `results/gate2/`, so the rest of the pipeline runs without an API key
- Node.js 18+ with the `docx` package — only for `scripts/05_report/25_manuscript_docx.js`
- About 20 GB of disk for the recount3 cache; the full reanalysis took roughly 6 h on 12 cores, the induced-imbalance experiment and the repeated permutations a further day on 40 cores

Set the project path at the top of the shell scripts (`PROJ=`); the R scripts read the same constant or the `COV_AUDIT_PROJ` environment variable.

## Run order

| Step | Scripts | Output |
|---|---|---|
| 1 | `00_pilot/01*_download_recount3_metadata*.sh` → `02_parse_metadata_blood_candidates.R`, `02m_parse_mouse_metadata.R`, `02c_select_all_tissue_candidates.R` | project and sample metadata, candidate lists |
| 2 | `01_survey/02_build_design_packets.R` → `03_classify_design_llm.py` (API) → `04_eval_design_llm.py`, `05_make_review_sheet.py` | design verdicts, validation |
| 3 | `00_pilot/03m_extract_mouse_sex_genes.sh` → `04m_mouse_sex_inference.R` | mouse sex from expression |
| 4 | `02_reanalysis/22_make_release_config.R` | `config/release_{human,mouse}_runs.tsv` |
| 5 | `02_reanalysis/10_run_main_de.sh`, `10m_run_main_de_mouse.sh` (drivers in `00_pilot/03*_pilot_de*.R`) → `11b_collect_main_summary_v2.R`, `11c_collect_extra_models.R` | per-study DE tables; `results/summary/main_*_summary_v2.tsv` |
| 6 | `02_reanalysis/12_c2_gsea.R` → `13_c2_summary.R` | pathway-level summaries |
| 7 | `03_qc/20_extract_accession_levels.sh` → `21_pseudoreplication_audit.R`; `15_prevalence_null.R`; `14_sens_sex_mismatch.R`; `22_batch_imbalance_survey.R`; `25_delta_regression.R` | QC tables |
| 8 | `04_gate2/16_gate2_fetch.py` → `17_gate2_llm.py` (API) → `18_gate2_summary.R` | publication audit |
| 9 | `06_revision/43_within_group_control.R`, `45_within_group_pathway.R` | within-group control by class and bin |
| 10 | `06_revision/42_make_synthetic_list.R` → `41_run_synthetic.sh` (driver `40_synthetic_imbalance.R`) → `44_synthetic_summary.R`, `50_biasvar_summary.R` | induced-imbalance experiment |
| 11 | `06_revision/47_run_content.sh` (driver `46_content_of_change.R`) → `49_content_summary.R`, `51_v4_numbers.R` | repeated permutations, content of the change |
| 12 | `06_revision/48_tissue_heterogeneity.R`; `52_gtex_enrichment.R` → `53_gtex_summary.R`; `30_single_sex_stratum.R` … `39_mechanism_summary.R` | heterogeneity, GTEx enrichment, revision analyses |
| 13 | `05_report/26_draft_figures_v3.R`, `31_supp_figures_v3.R`, `32_figures_v4.R`, `27_source_data_xlsx.py` | figures and source data |

Permutation seeds are derived from the study accession, so every model of a study is reproducible from its accession alone.

## Data deposited separately

Per-study differential-expression tables for every configured study (`de_M0`, `de_M1`, `de_M0perm`, `de_M0permw` per study, gzip-compressed TSV, about 6 GB) are archived on Zenodo: [DOI to be added]. Each folder also holds the study's phenotype table (`pheno.tsv`) and summary row.

## Citation and license

Code is released under the MIT License; tables under `results/` and `config/` under CC BY 4.0. Please cite the paper and the archived version of this repository (see `CITATION.cff`).
