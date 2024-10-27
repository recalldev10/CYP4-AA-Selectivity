# CYP4 Implementation Guide

## Environment Setup

### Prerequisites
1. WSL (Windows Subsystem for Linux)
2. VPN access
3. Cluster access
4. [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
5. PyMOL
6. AMBER
7. Required Python packages:
   - RDKit
   - NumPy
   - Pandas

### Initial Setup Steps
1. Install WSL following Windows documentation
2. Set up VPN access
3. Configure cluster access
4. Install LocalColabFold following GitHub instructions

## Project Implementation

### 1. Sequence Preparation
1. Save sequence FASTA from UniProt as `[cyp].fasta.txt`
2. Modify ColabFold script to utilize single GPU
3. Generate predictions using AlphaFold:
```bash
for i in 16 32 64 128 256; do colabfold_batch --num-seeds 5 --nu
```

### 2. Structure Analysis
1. Open generated PDB files in PyMOL
2. Align to reference PDB (rec.pdb)
3. Upload heme and save as new structure
4. Dock AA.sdf to each structure
5. Analyze AA binding to determine optimal CYP binding

### 3. Data Analysis with RDKit
```python
import glob
import gzip
from rdkit import Chem
import numpy as np
import pandas as pd
from rdkit.Geometry import rdGeometry

# Set heme iron position
FEpos = rdGeometry.Point3D(63.62, 22.366, 56.693)

def getC(pose):
    for atom in pose.GetAtoms():
        if int(atom.GetDegree()) == 1 and atom.GetSymbol() == 'C':
            return atom

# Analysis script for processing poses
ascores = []
for fname in glob.glob('aa256/*.gz'):
    for i, mol in enumerate(Chem.ForwardSDMolSupplier(gzip.open(fname))):
        props = mol.GetPropsAsDict()
        C = getC(mol)
        Cpos = mol.GetConformer().GetAtomPosition(C.GetIdx())
        dist = np.linalg.norm(Cpos-FEpos)
        props['Cdist'] = dist
        props['fname'] = fname
        props['idx'] = i
        ascores.append(props)

# Create and filter DataFrame
df = pd.DataFrame(ascores)
filtered_df = df[df['Cdist'] < 4]
sorted_df = filtered_df.sort_values(by=['CNNaffinity'], ascending=True)
top_three_poses = sorted_df.head(3)
```

### 4. Molecular Simulations

For each PDB file:

1. **Directory Setup**
   - Create separate directory for each complex
   - Unpack `files.tgz` in same directory as PDB

2. **Complex Preparation**
   - Run `sed` commands to combine AA and CYP-heme complex
   - Separate AA ligand from complex in PyMOL:
     - Select AA → Extract object
     - Save AA as `aa.pdb`
     - Save CYP-heme as `complex.pdb`

3. **Cluster Processing**
   - Transfer directory to cluster
   - Install `preamber.py` from Koes AMBER GitHub
   - Load AMBER module
   - Initialize GPU node:
     ```bash
     srun --gres=gpu:1 --pty -p dept_gpu /bin/bash -i
     ```

4. **AMBER Setup**
   - Run preamber.py: `/path/to/preamber.py -s amber.py`
   - Combine `aa_amber.pdb` with `complex.pdb`
   - Remove CONNECT lines:
     ```bash
     sed -ibkup '/CONECT/d' complex.pdb
     ```
   - Modify `tleap.in` for bond errors (configuration used for AMBER)
   - Run `tleap.in`

5. **Equilibration and Analysis**
   - Run equilibration commands
   - Convert trajectory:
     ```bash
     cpptraj -p complex.prmtop
     trajin complex_md2.nc
     trajout converted_trajectory.dcd dcd
     run
     quit
     ```
   - Analyze in PyMOL:
     ```bash
     pymol complex.prmtop complex_md2.dcd
     ```

## Notes
- Verify GLG and HEM residue connections before production MDs
- Consult with Dr. Koes before proceeding to production MDs
- Additional production MD steps will be provided when reaching that stage

## Troubleshooting
- For bond errors, check and modify `tleap.in`
- Verify residue numbering when combining structures
- Monitor GPU allocation on cluster
