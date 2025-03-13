# 🧬 Computational Structural Biology Toolkit

Welcome to my collection of Python scripts for analyzing and processing protein structures! Whether you're working with **PDB files**, calculating **RMSD**, or analyzing **pi-stacking interactions**, this repository has you covered. 🚀

---
## 📂 **Scripts Overview**

### 1️⃣ **Process_PDB.py** 🏗️
**What it does:**
- Processes PDB files to extract relevant structural information.
- Can filter, modify, and organize atomic data.

**How to use:**
```bash
python Process_PDB.py <input_pdb> <output_pdb> [options]
```

---
### 2️⃣ **Pi_Stacking_Analysis.py** 📏
**What it does:**
- Computes **pi-stacking interactions** between aromatic residues in a given set of PDB files.
- Outputs **distances** and **angles** between ring systems.

**How to use:**
```bash
python Pi_Stacking_Analysis.py <pdb_directory> <output_csv> <atom_list>
```

---
### 3️⃣ **Split_PDBs_Ensemble.py** ✂️
**What it does:**
- Extracts individual models from **ensemble PDB files**.
- Saves each model as a separate PDB file.

**How to use:**
```bash
python Split_PDBs_Ensemble.py <input_pdb>
```

---
### 4️⃣ **RMSD.py** 🔄
**What it does:**
- Computes **Root Mean Square Deviation (RMSD)** between two structures.
- Supports **backbone-only** or **all-atom RMSD**.

**How to use:**
```bash
python RMSD.py <structure_A.pdb> <structure_B.pdb> [--backbone]
```

---
### 5️⃣ **Calculate_Deviation_Diversity_Ens.py** 📊
**What it does:**
- Computes **ensemble deviation** (RMSD vs. reference crystal structure).
- Computes **ensemble diversity** (pairwise RMSD between ensemble members).

**How to use:**
```bash
python Calculate_Deviation_Diversity_Ens.py <crystal_structure.pdb> <ensemble_dir> <output_file.csv>
```

---
### 6️⃣ **Analyse_AlphaFold_Outputs.py** 🤖
**What it does:**
- Analyzes **AlphaFold-generated PDBs**.
- Extracts **pLDDT scores** and computes **backbone RMSD**.

**How to use:**
```bash
python Analyse_AlphaFold_Outputs.py <input_pdb_dir> <output_csv> [ref_pdb]
```
```

---
<a href="https://www.linkedin.com/in/niayesh-zarifi/" target="_blank">
  <img height="30" style="vertical-align: middle;" src="https://cdn4.iconfinder.com/data/icons/colorful-guache-social-media-logos-1/159/social-media_linkedin-512.png" alt="LinkedIn"/>
</a>

---
✨ *Happy coding!* ✨

