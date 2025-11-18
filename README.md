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

Your working directory should look like this:

```
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
```

- The script automatically detects the `models/` folder as long as it is in the **same directory** as the script.
- The **GFF3** and **parent IDs file** must be provided as command-line arguments.

---

## 🧩 Input Files

### 1. GFF3 file
Contains genomic annotations.

Example:  
`merged_gff_for_proteomics.gff3`

### 2. Parent Gene IDs
Text file with one parent gene ID per line.  
These must match the identifiers in your GFF.

Example `parent_gene_ids.txt`:

```
LOC_Os01g01484
Os11t0537400
OsNip_01G000390
```

### 3. Structural Models
All PDB/CIF files must be stored under:

```
models/
```

The script automatically finds files whose names match the gene IDs in the GFF.

---

## ▶️ Running the Script

From the folder containing `eshan_v8_cbf.py`, run:

```
python eshan_v8_cbf.py merged_gff_for_proteomics.gff3 parent_gene_ids.txt
```

- Uses all available CPU cores by default
- Prints debug information such as:
  - models directory path
  - number of CPU cores used
  - which models were matched or missing

---

## 📝 Output

The script produces:

- A consolidated annotation file (CSV/TSV depending on script version)
- Log messages for each processed gene
- Warnings for:
  - missing gene models  
  - GFF entries without parents  
  - files that do not match expected gene IDs  

The output file typically includes:

- gene ID  
- transcript ID  
- residue index  
- amino acid  
- structure file used  
- coordinates  
- additional annotations

---

## 🧪 Testing

To test with only a few genes:

```
head -n 3 parent_gene_ids.txt > test_ids.txt
python eshan_v8_cbf.py merged_gff_for_proteomics.gff3 test_ids.txt
```

---

## 🛠 Installation & Requirements

Requires Python **3.8+** and the following packages:

```
pip install biopython pandas tqdm
```

Optional (if needed):

```
pip install numpy
```

---

## 🐛 Troubleshooting

### Model marked missing, but file exists?
Check:

- The parent ID matches the filename prefix **exactly**
- No extra whitespace in `parent_gene_ids.txt`
- File extensions are `.pdb` or `.cif`

### All models missing (all red)?
Most common cause:
- Wrong parent IDs file  
- IDs don’t match the GFF  
- Hidden carriage returns (`\r`) in the text file  

### Multiprocessing crash?
Try forcing single-core mode:

```
export OMP_NUM_THREADS=1
python eshan_v8_cbf.py merged_gff_for_proteomics.gff3 parent_gene_ids.txt
```


## 📬 Support

If you encounter issues, feel free to open an Issue or contact the contributor who provided the script.
