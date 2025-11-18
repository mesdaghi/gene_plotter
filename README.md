# gene_plotter
Finds matching gene models in the GFF Matches them to structural files in ./models/ Extracts sequences, coordinates, residues, etc. Generates a combined output table Runs using multiple CPU cores (auto-detected)
# Eshan Structural Annotation Pipeline

This repository contains the script **`eshan_v8_cbf.py`**, which processes rice (Oryza sativa) gene models using GFF3 annotations and corresponding structural models (PDB/CIF).  
It maps GFF features to protein structures, extracts residue-level information, and generates a unified annotation table.

---

## ✨ Features

- Parses **GFF3** gene models
- Maps each gene ID to its **AlphaFold / structural model** (`.pdb` or `.cif`)
- Supports multiple model naming conventions (LOC_Os..., OsNip..., OsXt..., Afu..., etc.)
- Handles splice isoforms and versioned files
- Extracts:
  - sequences  
  - coordinates  
  - residue-level annotations  
- Uses **multiprocessing** for speed
- Produces a clean merged output table

---

## 📁 Directory Structure

Your working directory should look like:

project/
├── eshan_v8_cbf.py
├── merged_gff_for_proteomics.gff3
├── parent_gene_ids.txt
└── models/
├── LOC_Os01g01484.1.pdb
├── LOC_Os01g01484.4.pdb
├── Os11t0591000-01.pdb
├── OsNip_01G000390_03.pdb
└── ... more structures ...


- The script automatically detects the `models/` folder as long as it is in the **same directory** as the script.
- The **GFF3** and **parent IDs file** must be provided as command-line arguments.

---

## 🧩 Input Files

### **1. GFF3 file**
Contains genomic annotations.

Example:  
`merged_gff_for_proteomics.gff3`

### **2. Parent Gene IDs**
Text file with one parent gene ID per line.  
These must match the identifiers in your GFF.

Example `parent_gene_ids.txt`:

LOC_Os01g01484
Os11t0537400
OsNip_01G000390


### **3. Structural Models**
All PDB/CIF files must be stored under:

./models


The script automatically finds files whose names match the gene IDs in the GFF.

---

## ▶️ Running the Script

From the folder containing `eshan_v8_cbf.py`, run:

```bash
python eshan_v8_cbf.py merged_gff_for_proteomics.gff3 parent_gene_ids.txt

