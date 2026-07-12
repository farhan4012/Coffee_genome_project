# Coffee Genome Project

An end-to-end genomics pipeline for variant identification and candidate gene annotation in allotetraploid *Coffea arabica*, built around the ET-39 reference genome (CC and CE subgenomes).

## Traits Covered
- Low Caffeine
- Coffee Leaf Rust Resistance (SH3)
- Coffee Berry Disease Resistance (Ck-1)
- Dwarfism / Compactness
- Giant Bean Size

## Data Setup
GFF3 files are not included due to file size. Place them in `data/gff3/` before running and include the path in app.py file at line 40-41.

## Run
pip install -r requirements.txt
python app.py

Open browser at http://localhost:5000

## Stack
Flask, minimap2, bcftools, Chart.js

## Author
Md Farhan Alam, B.Tech Biotechnology, IARI Pusa 2026
