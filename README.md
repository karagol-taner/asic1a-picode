# Position and charge change predict measured pH50 shifts in a pI-code ladder test of acid-sensing ion channel 1a

---

## Abstract

Protein design methods that restrict substitutions to single base changes assume that the resulting near-isosteric swaps are the ones least likely to disturb a fold, an assumption supported so far by side-chain geometry and by substitution costs drawn from variant-effect predictors and curated human variation. None of those is a measurement of function. We test it against acid-sensing ion channel 1a, a target carrying a measured, continuous, per-substitution functional readout, using 117 single substitutions over 38 positions recorded under one protocol in one laboratory, with an independent series of 33 substitutions supplying two ladder edges the first set cannot reach. What predicts the measured shift is a property of the position and not of the residue pair. Relative side-chain burial orders it at -0.595 across the 38 scanned positions (permutation p = 0.0001) on the acidified structure, falling to -0.111 (p = 0.508) on the resting one, so the descriptor is conformer-specific. Charge change carries the measured effect at +0.170 pH units per charge unit, and over the eight scorable edges charge-change magnitude orders the shifts at +0.914 (p = 0.004) while substitution cost orders them in the predicted direction without reaching significance (Spearman rho -0.667, exact permutation p = 0.083). Every acidic-pocket carboxylate pair is intra-subunit across four assemblies spanning pH 5.7 to 8.5 while the cavity they line is completed by the neighbouring subunit. Holding position, charge change and side-chain volume fixed across 101 substitutions at 22 positions, ladder membership changes the activation shift by -0.035 pH units (permutation p = 0.437) and the steady-state desensitization shift by -0.001 units (p = 0.974). A sign test over charge-matched pairs at single positions does not separate them either (p = 0.110), and reaches nominal significance only when the indicator is widened past the twelve edges the method defines. Codon-level constraint therefore delivers charge bookkeeping with no detectable tolerance effect, and substitution costs derived from pathogenicity should not be quoted as predicting pH behaviour.

Keywords: pH-responsive protein design, Substitution tolerance, Acid-sensing ion channel, Variant effect prediction, Codon-level constraint, Proton sensing

---

## Repository layout

```
├── README.md                        This file
├── codes/
│   ├── picode_paths.py              Resolves the ladder and burial constants
│   ├── fetch_structures.py          Downloads the four assemblies from the PDB
│   ├── asic_ph50_dataset.py         Primary cohort from the publisher PDF, cross-checked
│   ├── asic_ph50_canessa.py         Second cohort from Europe PMC full-text XML
│   ├── asic_numbering_map.py        Chicken, human and rat numbering by alignment
│   ├── asic_edge_coverage.py        Which ladder edges the measured set covers
│   ├── asic_fixed_effects.py        Primary test: position fixed effects, charge, volume
│   ├── asic_second_shell.py         Electrostatic determinants and their calibration
│   ├── asic_results.py              Every reported number, computed from data/ and structures/
│   ├── extract_edge_costs.py        Per-edge costs; needs the unpublished scoring code
│   ├── extract_interface_carboxylates.py  Buried uncompensated carboxylates, both frames
│   ├── asic_edge_vs_cost.py         Eight-edge rank test; needs the unpublished scoring code
│   ├── export_ladder_definition.py  Refreshes data/*.json from the scoring code
│   ├── build_figures.py             Figures 1-6 and S1-S2
│   ├── build_md_figure.py           Supplementary Figure S3 from the simulation reports
│   ├── make_md_notebook.py          Regenerates the notebooks in md/
│   ├── check_notebooks.py           Every name resolves on a Run All
│   └── test_md_notebooks.py         Notebooks against their templates
├── data/
│   ├── asic1a_ph50.csv              117 substitutions plus wild type, both readouts, SEM, n
│   ├── asic1a_ph50_canessa.csv      33 substitutions, second laboratory, never pooled
│   ├── ladder_definition.json       Ladder edges and per-residue charges
│   ├── burial_thresholds.json       Accessibility cutoffs and the Tien maximum-ASA table
│   └── sources/                     Publisher supplementary files, full-text XML, UniProt
├── results/
│   ├── edge_costs.csv               Per-edge substitution costs
│   └── interface_carboxylates.csv   Monomer and assembly frames per deposited assembly
├── structures/                      PROPKA output for the low-pH and resting assemblies
├── md/                              MD notebooks, CHARMM-GUI inputs, GROMACS parameters,
│                                    and the per-frame analysis behind Figure S3
├── figures/                         Figures 1-6 and S1-S3, PNG at 400 dpi
└── supplementary/                   Supplementary Information and the data workbook
```

---

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `data/asic1a_ph50.csv` | 118 | Primary cohort: 117 single substitutions over 38 positions plus the wild type, activation and steady-state desensitization pH50 with SEM and n, structural domain. |
| `data/asic1a_ph50_canessa.csv` | 34 | Second cohort: 33 substitutions over 4 positions from an independent laboratory, plus the wild type. Never pooled with the first. |
| `results/edge_costs.csv` | 25 | Per-edge substitution cost for the twelve core, five arginine and eight proposed-extension edges, by Cliff's delta and by pathogenic fraction. |
| `results/interface_carboxylates.csv` | 4 | Buried charge-uncompensated carboxylates per deposited assembly, scored in the monomer frame and in the assembly frame. |

Both cohorts are parsed from the publisher files in `data/sources/` at run time rather than transcribed. The primary cohort is Supplementary Table 2 of Bignucolo et al. (2025), cross-checked against the four-decimal fitted values in the same study's Supplementary Data 1; the second is Table 1 of Wu et al. (2019).

---

## Reproducing the analysis

```bash
pip install numpy scipy matplotlib biopython openpyxl pdfplumber

python3 codes/fetch_structures.py            # the four assemblies, from the PDB
python3 codes/picode_paths.py                # where the constants come from
python3 codes/asic_ph50_dataset.py data/sources/42003_2025_9090_MOESM2_ESM.pdf \
    --xlsx data/sources/42003_2025_9090_MOESM4_ESM.xlsx -o data/asic1a_ph50.csv
python3 codes/asic_ph50_canessa.py data/sources/PMC6497441.xml \
    -o data/asic1a_ph50_canessa.csv
python3 codes/asic_numbering_map.py \
    --chick data/sources/up_Q1XA76.fasta --human data/sources/up_P78348.fasta
python3 codes/asic_fixed_effects.py data/asic1a_ph50.csv \
    --props data/sources/42003_2025_9090_MOESM4_ESM.xlsx --readout act
python3 codes/asic_edge_coverage.py data/asic1a_ph50.csv
python3 codes/asic_second_shell.py --pka structures/9E4F.pka \
    --pka-resting structures/9E4A.pka --ph50 data/asic1a_ph50.csv
python3 codes/asic_results.py                # every reported number
python3 codes/build_figures.py               # Figures 1-6 and S1-S2
```

Every script runs on the files in this repository once the assemblies are fetched, and exits non-zero on a fatal condition rather than continuing. `codes/asic_fixed_effects.py --edges implemented` refits the primary model on the wider edge set reported as a sensitivity.

---

## Software dependencies

| Domain | Packages |
| --- | --- |
| Data wrangling | `numpy` |
| Plotting | `matplotlib` |
| Statistics | `scipy` |
| Structural / sequence | `biopython` (`Bio.PDB`, `Bio.Align`) |
| File formats | `openpyxl`, `pdfplumber` |
| Molecular dynamics | GROMACS and CHARMM-GUI, used through the notebooks in `md/` |

Tested on Python 3.10.

---

## Data availability

* The four biological assemblies, 2QTS, 6VTL, 9E4A and 9E4F, are public depositions and are downloaded from the Protein Data Bank by `codes/fetch_structures.py`, which checks each one against the chain count, atom count and residue span the analyses were run against.
* Measurement data are third-party and published: Bignucolo et al. (2025), *Communications Biology* 8:1692, and Wu et al. (2019), *eLife* 8:e45851. The source files they were parsed from are included under `data/sources/`.
* UniProt records are from release 2026_03 (`https://rest.uniprot.org/`); the maximum-accessibility reference is that of Tien et al. (2013), *PLoS ONE* 8:e80635.
* The scoring implementation of the design method under test is the subject of a patent application and is not released. `extract_edge_costs.py`, `extract_interface_carboxylates.py` and `asic_edge_vs_cost.py` execute it and so do not run from a clone; all three exit non-zero saying so, and the output of the first two is committed under `results/`. The constants the rest of the analysis needs are committed in `data/` and are published in Figure 1, Table 2 and the Methods of the manuscript.

---

## Citation

If you use this code or dataset, please cite:

> Karagöl T., Karagöl A. *Position and charge change predict measured pH50 shifts in a pI-code ladder test of acid-sensing ion channel 1a.* (2026).

A BibTeX entry will be added upon publication of the journal record.

---

## Contact

Questions and issues are welcome, please open a GitHub Issue or contact the corresponding authors.

---

## License

The code in this repository is released under the MIT License. Curated datasets and figures are released under CC-BY 4.0. The underlying measurement data, Protein Data Bank entries and UniProt records are subject to the licenses of those upstream providers.
