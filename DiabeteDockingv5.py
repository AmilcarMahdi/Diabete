import vcf, os, subprocess, warnings, re
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from pocket_finder import trouver_poche_profonde 

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

# --- FONCTION DE NETTOYAGE CRITIQUE ---
def reparer_types_atomes(pdbqt_file):
    """Supprime les symboles de charge comme '+' ou '-' dans la colonne des types d'atomes"""
    with open(pdbqt_file, "r") as f:
        lines = f.readlines()
    
    with open(pdbqt_file, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                # La colonne du type d'atome est à la fin (77-79)
                type_part = line[77:80]
                if "+" in type_part or "-" in type_part:
                    # On remplace N1+ par N, C1+ par C, etc.
                    nouveau_type = type_part.replace("+", "").replace("-", "").replace("1", "")
                    line = line[:77] + nouveau_type.ljust(3) + "\n"
            f.write(line)

def preparer_ligand_pro(name, smiles):
    out = f"{name}.pdbqt"
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        prepper = MoleculePreparation()
        prepper.prepare(mol)
        with open(out, "w") as f: f.write(prepper.write_pdbqt_string())
        
        # --- RÉPARATION ICI ---
        reparer_types_atomes(out)
        return out
    except: return None

def criblage_vina(rsid, info):
    meilleur_score = 0.0
    meilleure_mol = "Aucun"
    center = trouver_poche_profonde(info['pdb'])
    cx, cy, cz = center

    for mol_nom, smiles in LIGAND_LIBRARY.get(info['gene'], []):
        ligand_file = preparer_ligand_pro(mol_nom, smiles)
        conf = f"conf_{rsid}.txt"
        with open(conf, "w") as f:
            size = 14.0 if info['pdb'] == "6C3O" else 18.0
            f.write(f"receptor = {info['pdb']}.pdbqt\nligand = {ligand_file}\n")
            f.write(f"center_x={cx:.3f}\ncenter_y={cy:.3f}\ncenter_z={cz:.3f}\n")
            f.write(f"size_x={size}\nsize_y={size}\nsize_z={size}\nexhaustiveness=8\n")

        print(f"      - Calcul Vina : {mol_nom}...")
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        
        scores = re.findall(r"-\d+\.\d+", res.stdout)
        if scores:
            score = float(scores[0])
            print(f"        -> ✅ Score : {score} kcal/mol")
            if score < meilleur_score:
                meilleur_score = score
                meilleure_mol = mol_nom
        elif res.stderr:
            print(f"        ⚠️ ERREUR : {res.stderr.strip()[:60]}")
            
    return meilleure_mol, meilleur_score

# --- (Garde le reste du script main et les dictionnaires identiques) ---
