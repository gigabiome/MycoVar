# 🧬 MycoVar

**MycoVar** is a toolkit for performing **phylogenetics**, **lineage marker analysis**, and **resistotyping** on *Mycobacterium tuberculosis* (Mtb).  
It provides a streamlined workflow through Jupyter notebooks, backed by a shared functions library and curated reference databases.

---

## 📂 Repository Structure

MycoVar/
├── phylo.ipynb # Phylogenetic analysis
├── lineage.ipynb # Lineage marker analysis
├── resisto.ipynb # Resistotyping analysis
├── functions.py # Core helper functions used across notebooks
├── bin/ # Databases and reference files
├── input/ # Place your input data here
├── output/ # Results will be saved here

yaml
Copy
Edit

---

## ⚡ Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/your-username/MycoVar.git
cd MycoVar
2. Install dependencies
It is recommended to use a fresh virtual environment.
Required packages include Python 3.x, Jupyter, and common bioinformatics/data analysis libraries.

bash
Copy
Edit
# Example using pip
pip install -r requirements.txt
(If you don’t yet have a requirements.txt, you can generate one with pip freeze > requirements.txt.)

3. Run notebooks
Launch Jupyter and open any of the analysis notebooks:

bash
Copy
Edit
jupyter notebook
phylo.ipynb → Build phylogenetic trees

lineage.ipynb → Identify lineage markers

resisto.ipynb → Detect resistance profiles

All results are automatically saved to the output/ folder.

🧩 Workflow Overview
Place input files (e.g. variant call data or sequence data) into input/.

Open the relevant notebook (phylo, lineage, or resisto).

Execute the cells step-by-step.
Each notebook calls functions from functions.py to standardize analysis.

Results are written to output/.

🤝 Contributing
Contributions, improvements, and suggestions are welcome!
To contribute:

Fork the repo

Create a feature branch (git checkout -b feature-name)

Commit and push changes

Open a Pull Request

📜 License
This project is licensed under the MIT License.
Feel free to use and adapt with attribution.

🙏 Acknowledgements
Built for research on Mycobacterium tuberculosis phylogeny, lineage, and resistance.

Uses curated databases included in the bin/ directory.

