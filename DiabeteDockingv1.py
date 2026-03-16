import vcf, os, subprocess, warnings
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

TABLEAU_STRATEGIQUE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'box': [11.0, 14.5, 19.0], 'drug': 'Insulin_Mimetic'},
    'rs1799884': {'gene': 'GCK', 'pdb': '1V4S', 'box': [18.5, 10.0, 12.5], 'drug': 'GCK_Activator'},
    'rs11212617': {'gene': 'ATM', 'pdb': '7SAY', 'box': [162.0, 175.0, 160.0], 'drug': 'Metformine'},
    'rs5219': {'gene': 'KCNJ11', 'pdb': '6C3O', 'box': [115.0, 110.0, 95.0], 'drug': 'Gliclazide'},
    'rs4124874': {'gene': 'SLC22A1', 'pdb': '8H66', 'box': [125.0, 128.0, 142.0], 'drug': 'Bypass_Molecule'},
    'rs7903146': {'gene': 'TCF7L2', 'pdb': '2G4B', 'box': [25.0, 15.0, -8.0], 'drug': 'GLP1_Analogue'}
}

MOLECULES_SMILES = {
    'Metformine': 'CN(C)C(=N)N=C(N)N',
    'Gliclazide': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3',
    'Insulin_Mimetic': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O',
    'GCK_Activator': 'CC1=CC=C(C=C1)S(=O)(=O)N',
    'Bypass_Molecule': 'CN1C2=C(C(=O)N(C1=O)C)N=CN2',
    'GLP1_Analogue': 'C1=CC=C(C=C1)CC(C(=O)O)N'
}

def nettoyage_vina_pro(file_path):
    """Nettoyage compatible Vina 1.2 : préserve ATOM et ajoute TORSDOF si besoin"""
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    clean = []
    has_torsdof = False
    for l in lines:
        if l.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
            clean.append(l)
            if "TORSDOF" in l: has_torsdof = True
            
    if not has_torsdof and "ROOT" in "".join(clean):
        clean.append("TORSDOF 0\n")
        
    with open(file_path, "w") as f:
        f.writelines(clean)

def preparer_ligand_v5_7(name, smiles):
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

def lancer_docking(rsid, info):
    ligand = preparer_ligand_v5_7(info['drug'], MOLECULES_SMILES[info['drug']])
    receptor = f"{info['pdb']}.pdbqt"
    if not os.path.exists(receptor): return "PDBQT manquant"
    
    # On nettoie le récepteur une fois pour toutes
    nettoyage_vina_pro(receptor)

    cx, cy, cz = info['box']
    conf = f"conf_{rsid}.txt"
    with open(conf, "w") as f:
        f.write(f"receptor = {receptor}\nligand = {ligand}\n")
        f.write(f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n")
        f.write("size_x = 22.0\nsize_y = 22.0\nsize_z = 22.0\nexhaustiveness = 8\n")

    try:
        # On lance Vina et on capture le flux de sortie
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        if res.returncode != 0: return f"Vina Error"
        
        # Extraction du score réel (ex: -7.4)
        for line in res.stdout.split('\n'):
            if "   1 " in line:
                parts = line.split()
                return parts[1] # On prend la 2ème colonne (le score)
    except: return "Erreur"
    return "N/A"

def main():
    print("="*90)
    print("  DIABETE-GENE V5.7 : EXTRACTION DU SCORE RÉEL & FIX BOX")
    print("="*90)
    
    import vcf
    if not os.path.exists("patient_expert.vcf"): return
    reader = vcf.Reader(filename="patient_expert.vcf")
    for record in reader:
        if record.ID in TABLEAU_STRATEGIQUE:
            info = TABLEAU_STRATEGIQUE[record.ID]
            print(f"\nMUTATION: {record.ID} ({info['gene']})")
            score = lancer_docking(record.ID, info)
            status = "✅ SUCCÈS" if "-" in str(score) else "⚠️ ÉCHEC (Vérifier Box)"
            print(f"   AFFINITÉ : {score} kcal/mol | {status}")

if __name__ == "__main__":
    main()
