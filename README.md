# Protein-Mutation-Feature-Extractor-Work-in-Progress-
This project aims to extract structured features for protein residues to train machine learning models that predict the functional impact of protein mutations.
Features Implemented So Far
	•	Fetch data from PDB and ClinVar.
	•	Extract structural and functional features per residue:
	◦	Solvent accessibility
	◦	Clinical significance
	◦	Mutation types (protein-altering: missense, frameshift)
	◦	Secondary structure
	•	Sequence conservation calculation using UniProt API and multiple sequence alignment (MAFFT via Biopython).

Planned Features
	•	Streamlined pipeline combining all modules.
	•	Automated variant-to-structure mapping.
	•	Full ML-ready dataset generation for mutation effect prediction.
	•	Additional conservation scoring improvements.

Project Structure

├── conservation_calculator  # Compute conservation scores
├── data_fetch               # Modules for fetching PDB, ClinVar, UniProt data
├── results                  # Temporary files and alignment outputs
├── structural_features      # Extract structural features from PDB files
├── variant_mapping          # Map variants to protein residues (planned)
├── requirements             # Python dependencies

![Work in Progress](https://img.shields.io/badge/status-wip-orange)
