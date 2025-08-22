# 🧬 MycoVar

**MycoVar** is a toolkit for performing **phylogenetics**, **lineage marker analysis**, and **resistotyping** on *Mycobacterium tuberculosis* (Mtb).  
It provides a streamlined, notebook-driven workflow backed by a shared function library and bundled reference databases.

---

## ✨ Features

- **Phylogenetics**: build trees and explore evolutionary relationships.  
- **Lineage analysis**: detect lineage markers from variants.  
- **Resistotyping**: infer resistance profiles from genomic data.  
- **Notebook-first**: three focused Jupyter notebooks that call a common `functions.py`.  
- **Self-contained**: curated databases stored under `bin/`; inputs in `input/`, results in `output/`.

---

## 📂 Repository Structure

```
MycoVar/
├── phylo.ipynb           # Phylogenetic analysis
├── lineage.ipynb         # Lineage marker analysis
├── resisto.ipynb         # Resistotyping analysis
├── functions.py          # Shared helper functions used by all notebooks
├── bin/                  # Databases and reference files (read-only)
├── input/                # Place your input data here
└── output/               # All results are written here
```

> 💡 Run notebooks from the repository root so they can import `functions.py` and find `bin/`.

---

## 🚀 Quick Start

### 1) Clone the repository
```bash
git clone https://github.com/your-username/MycoVar.git
cd MycoVar
```

### 2) Launch Jupyter and run an analysis
```bash
jupyter notebook
```

Open one of:
- **`phylo.ipynb`** → build/visualize phylogenetic trees  
- **`lineage.ipynb`** → call lineage markers  
- **`resisto.ipynb`** → infer resistance profiles  

All outputs are written to the **`output/`** directory.

---

## 📥 Inputs

Place the required input files in `input/`.

---

## 📤 Outputs

Each notebook writes its results to `output/`, typically including:
- **Phylo**: tree files, summary tables, and optional figures  
- **Lineage**: per-sample lineage calls  
- **Resisto**: per-sample resistotype predictions and supporting evidence tables

---

## 🔧 How It Works

- The notebooks import shared utilities from **`functions.py`** to standardize data loading, QC, marker lookup, and result writing.  
- Reference databases (e.g., lineage markers, resistance loci) are read from **`bin/`**. Avoid renaming files in `bin/` unless you also update the corresponding paths in `functions.py` or the notebooks.

---

## 📦 Dependencies

Minimum recommended stack (unversioned; pin as needed):
- Python ≥ 3.9  
- Jupyter / jupyterlab  
- numpy, pandas, scipy  
- matplotlib (and/or plotly/altair for viz if used in notebooks)  
- biopython  
- scikit-learn (if used for clustering or models)  
- scikit-bio or ete3 (if used for tree handling)  

Install via:
```bash
pip install -r requirements.txt
```

**Example `requirements.txt` (adapt as needed):**
```
jupyter
pandas
numpy
scipy
matplotlib
biopython
scikit-learn
scikit-bio
ete3
```

---

## 🙏 Acknowledgements

- Developed for research on *Mycobacterium tuberculosis* phylogeny, lineage, and resistance at Centre For Tuberculosis NICD 
- Uses curated databases included in the `bin/` directory.  
- Thanks to the open-source bioinformatics community for foundational libraries.

---

## 📣 Citation

If you use **MycoVar** in your work, please cite this repository.  
(Add a DOI or citation block here once available.)
