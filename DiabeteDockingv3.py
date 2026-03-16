import vcf, os, subprocess, warnings
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
# Assure-toi que ton fichier s'appelle bien pocket_finder.py ou adapte le nom ici
from pocketfinderv12 import trouver_poche_profonde 

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

# --- BASE DE CONNAISSANCES MISE À JOUR (VRAIS MÉDICAMENTS) ---
TABLEAU_STRATEGIQUE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'drug': 'Insulin_Mimetic', 'smiles': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O'},
    'rs1799884': {'gene': 'GCK', 'pdb': '1V4S', 'drug': 'GCK_Activator', 'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N'},
    'rs11212617': {'gene': 'ATM', 'pdb': '7SAY', 'drug': 'Metformine', 'smiles': 'CN(C)C(=N)N=C(N)N'},
    'rs5219': {'gene': 'KCNJ11', 'pdb': '6C3O', 'drug': 'Gliclazide', 'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3'},
    'rs4124874': {'gene': 'SLC22A1', 'pdb': '8H66', 'drug': 'Bypass_Molecule', 'smiles': 'CN1C2=C(C(=O)N(C1=O)C)N=CN2'},
    'rs7903146': {'gene': 'TCF7L2', 'pdb': '2G4B', 'drug': 'GLP1_Analogue', 'smiles': 'C1=CC=C(C=C1)CC(C(=O)O)N'}
}

def nettoyage_vina_pro(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    clean = [l for l in lines if l.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"))]
    if not any("TORSDOF" in l for l in clean) and "ROOT" in "".join(clean):
        clean.append("TORSDOF 0\n")
    with open(file_path, "w") as f:
        f.writelines(clean)

def preparer_ligand_expert(name, smiles):
    out = f"{name}.pdbqt"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    prepper = MoleculePreparation()
    prepper.prepare(mol)
    with open(out, "w") as f:
        f.write(prepper.write_pdbqt_string())
    nettoyage_vina_pro(out)
    return out

def lancer_docking_v3(rsid, info):
    # 1. Préparer le VRAI médicament associé à la mutation
    ligand = preparer_ligand_expert(info['drug'], info['smiles'])
    
    receptor_pdbqt = f"{info['pdb']}.pdbqt"
    if not os.path.exists(receptor_pdbqt): return "PDBQT manquant"
    
    # 2. IA Pocket Finder
    coords = trouver_poche_profonde(info['pdb'])
    if not coords: return "Erreur PocketFinder"
    cx, cy, cz = coords

    # 3. Config Vina (Taille réduite pour plus de précision)
    conf = f"conf_{rsid}.txt"
    size = 18.0 if info['pdb'] != "6C3O" else 12.0 # On réduit pour les gros complexes
    with open(conf, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\nligand = {ligand}\n")
        f.write(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n")
        f.write(f"size_x = {size}\nsize_y = {size}\nsize_z = {size}\n")
        f.write("exhaustiveness = 12\n") # On augmente la puissance de calcul

    # 4. Exécution
    try:
        print(f"   [VINA] Docking Expert : {info['drug']} sur {info['pdb']}...")
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        for line in res.stdout.split('\n'):
            if "   1 " in line:
                return line.split()[1]
    except: return "Erreur"
    return "N/A"

def main():
    print("="*90)
    print("  DIABETE-GENE V3 : DOCKING DE PRÉCISION (VRAIS MÉDICAMENTS)")
    print("="*90)
    
    import vcf
    reader = vcf.Reader(filename="patient_expert.vcf")
    for record in reader:
        if record.ID in TABLEAU_STRATEGIQUE:
            info = TABLEAU_STRATEGIQUE[record.ID]
            print(f"\nANALYSE : {record.ID} ({info['gene']})")
            score = lancer_docking_v3(record.ID, info)
            print(f"   💊 MÉDICAMENT : {info['drug']} | AFFINITÉ : {score} kcal/mol")

if __name__ == "__main__":
    main()
