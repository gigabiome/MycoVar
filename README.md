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

### 2) Create an environment & install dependencies
Use your preferred environment manager. Example with **conda**:
```bash
conda create -n mycovar python=3.10 -y
conda activate mycovar
pip install -r requirements.txt
```

Or with **venv**:
```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

> If `requirements.txt` isn’t provided yet, see **Dependencies** below for a suggested list.

### 3) Launch Jupyter and run an analysis
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

Place the required input files in `input/`. Typical inputs may include:
- Variant calls (e.g., VCF), consensus FASTA, or tabular variant matrices
- Sample metadata (optional; CSV/TSV) for labeling/stratification

Refer to the first cell(s) of each notebook for the exact expected filenames/columns and any configurable parameters (paths, thresholds, marker sets).

---

## 📤 Outputs

Each notebook writes its results to `output/`, typically including:
- **Phylo**: tree files (e.g., Newick), summary tables, and optional figures  
- **Lineage**: per-sample lineage calls and confidence metrics  
- **Resisto**: per-sample resistotype predictions and supporting evidence tables

File names and formats are documented in the notebooks’ final cells.

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

> ⚠️ Only include libraries you actually use. If your notebooks rely on additional tools (e.g., plotly, seaborn), add them here.

---

## 🧪 Reproducibility

- Run notebooks from the repository root to ensure relative paths resolve.  
- Consider exporting executed notebooks to HTML for archival:
  ```bash
  jupyter nbconvert --to html --execute phylo.ipynb
  ```
- For strict reproducibility, pin package versions in `requirements.txt` and/or export a lock file:
  ```bash
  pip freeze > requirements.txt
  ```

---

## 🛠️ Troubleshooting

- **`ModuleNotFoundError: No module named 'functions'`**  
  Run Jupyter from the repo root (`MycoVar/`) or add the repo root to `PYTHONPATH`.

- **Databases not found**  
  Ensure `bin/` is present and paths in the notebook configuration cells match the filenames.

- **Permission errors writing outputs**  
  Verify you have write permissions to `output/` and that the directory exists.

---

## 🤝 Contributing

Contributions and suggestions are welcome!

1. Fork the repository  
2. Create a feature branch: `git checkout -b feature-name`  
3. Commit your changes: `git commit -m "Add feature"`  
4. Push to your fork and open a Pull Request

Please keep code style consistent and add/update notebook cells that document new parameters or outputs.

---

## 📜 License

This project is licensed under the **MIT License**. See `LICENSE` for details.

---

## 🙏 Acknowledgements

- Developed for research on *Mycobacterium tuberculosis* phylogeny, lineage, and resistance.  
- Uses curated databases included in the `bin/` directory.  
- Thanks to the open-source bioinformatics community for foundational libraries.

---

## 📣 Citation

If you use **MycoVar** in your work, please cite this repository.  
(Add a DOI or citation block here once available.)
