# HighHydrogenML
Python routines developed within the Marie Skłodowska-Curie Actions (MSCA) Postdoctoral Fellowship HighHydrogenML (grant agreement 101105610) for the publication: Martínez-Alonso, C., Vassilev-Galindo, V., Comer, B. M., Abild-Pedersen, F., Winther, K. T., & LLorca, J. (2024), Catalysis Science & Technology, 14(13), 3784-3799

**Get_Dataset.py**

Routine to create a pickle and a JSON file containing all relevant geometric, strain, and electronic descriptors (PSI, cell volume, weighted atomic radius, generalized coordination number, weighted electronegativity, weighted first ionization energy, outer electrons, and biaxial strain), as well as the adsorption energy directly from Quantum Espresso outputs. See Zenodo repository https://zenodo.org/doi/10.5281/zenodo.10364546 containing two ZIP folders (Bimetallic.zip and Pure_metals.zip) with Quantum Espresso files. The folders containing the Quantum Espresso files must have the exact same format as those within the ZIP folder for Get_dataset to work correctly. The routine uses GCN.py and PSI.py to compute the Generalized Coordination number and the Psi descriptors, respectively.

Usage:
_python Get_Dataset.py API.txt [Dataset.pickle]_

API.txt must contain the API key for using Materials Project within Python scripts (See https://next-gen.materialsproject.org/api for information)

Dataset.pickle is an already existing dataset file in which the user wants to append new entries for the dataset. Optional

**RandomForest.py** 

Train a Random Forest model for the prediction of adsorption energies. The dataset must be in xlsx format and the descriptor values must be already scaled using a MinMax scheme (see Zenodo repository https://zenodo.org/doi/10.5281/zenodo.10364546 for available datasets)

Usage:
_python RandomForest.py Dataset.xlsx_
