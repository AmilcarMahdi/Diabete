import vcf, os, subprocess, warnings
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
# Import de ton script de recherche de serrure
from pocket_finder import trouver_poche_profonde 

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

# --- BASE DE CONNAISSANCES STRATÉGIQUE (TRIAGE) ---
TABLEAU_STRATEGIQUE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'panne': 'Binding',    'calcul_vina': True,  'drug': 'Insulin_Mimetic', 'smiles': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O'},
    'rs5219':    {'gene': 'KCNJ11','pdb': '6C3O', 'panne': 'Binding',    'calcul_vina': True,  'drug': 'Gliclazide',      'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3'},
    'rs1799884': {'gene': 'GCK',   'pdb': '1V4S', 'panne': 'Activité',   'calcul_vina': False, 'drug': 'Activateur_GCK',  'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N'},
    'rs11212617':{'gene': 'ATM',   'pdb': '7SAY', 'panne': 'Activité',   'calcul_vina': False, 'drug': 'Metformine',      'smiles': 'CN(C)C(=N)N=C(N)N'},
    'rs7903146': {'gene': 'TCF7L2','pdb': '2G4B', 'panne': 'Expression', 'calcul_vina': False, 'drug': 'GLP1_Analogue',   'smiles': 'C1=CC=C(C=C1)CC(C(=O)O)N'},
    'rs4124874': {'gene': 'SLC22A1','pdb': '8H66', 'panne': 'Expression', 'calcul_vina': False, 'drug': 'Bypass_Molecule', 'smiles': 'CN1C2=C(C(=O)N(C1=O)C)N=CN2'}
}

def nettoyage_vina_pro(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    clean = [l for l in lines if l.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"))]
    if not any("TORSDOF" in l for l in clean) and "ROOT" in "".join(clean):
        clean.append("TORSDOF 0\n")
    with open(file_path, "w") as f:
        f.writelines(clean)

def preparer_ligand(name, smiles):
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

def lancer_docking_expert(rsid, info):
    # 1. Préparation du ligand
    ligand = preparer_ligand(info['drug'], info['smiles'])
    
    # 2. Recherche de la serrure via ton Pocket Finder
    print(f"   🔍 Analyse de la serrure profonde pour {info['pdb']}...")
    coords = trouver_poche_profonde(info['pdb'])
    if not coords: return "Erreur Géométrie"
    cx, cy, cz = coords

    # 3. Nettoyage du récepteur
    receptor_pdbqt = f"{info['pdb']}.pdbqt"
    if not os.path.exists(receptor_pdbqt): return "PDBQT manquant"
    nettoyage_vina_pro(receptor_pdbqt)

    # 4. Config et Exécution Vina
    conf = f"conf_{rsid}.txt"
    size = 18.0 if info['pdb'] != "6C3O" else 14.0
    with open(conf, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\nligand = {ligand}\n")
        f.write(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n")
        f.write(f"size_x = {size}\nsize_y = {size}\nsize_z = {size}\n")
        f.write("exhaustiveness = 12\n")

    try:
        print(f"   🚀 Lancement Vina Docking...")
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        for line in res.stdout.split('\n'):
            if "   1 " in line:
                return f"{line.split()[1]} kcal/mol"
    except: return "Erreur Système"
    return "N/A"

def main():
    print("="*95)
    print("  DIABETE-GENE V4 : SYSTÈME EXPERT DE DÉCISION (BINDING vs ACTIVITÉ vs EXPRESSION)")
    print("="*95)
    
    if not os.path.exists("patient_expert.vcf"):
        print("Erreur : patient_expert.vcf introuvable.")
        return
        
    reader = vcf.Reader(filename="patient_expert.vcf")
    for record in reader:
        if record.ID in TABLEAU_STRATEGIQUE:
            info = TABLEAU_STRATEGIQUE[record.ID]
            print(f"\n[DIAGNOSTIC] Mutation : {record.ID} | Gène : {info['gene']}")
            print(f"   Type de Panne : {info['panne'].upper()}")
            
            # --- LE FILTRE D'INNOVATION ---
            if info['calcul_vina']:
                print(f"   >>> Panne structurelle détectée. Simulation physique requise.")
                score = lancer_docking_expert(record.ID, info)
                print(f"   ✅ RÉSULTAT VINA : {score}")
            else:
                print(f"   >>> Panne non-structurelle. Calcul Vina ignoré (non pertinent).")
                print(f"   💊 STRATÉGIE CONSEILLÉE : {info['drug']}")

    print("\n" + "="*95)

if __name__ == "__main__":
    main()
