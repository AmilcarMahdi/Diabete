import vcf, os, subprocess, warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
# --- IMPORTATION DE TON SCRIPT POCKET FINDER ---
from pocketfinderv12 import trouver_poche_profonde 

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

TABLEAU_STRATEGIQUE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'drug': 'Insulin_Mimetic'},
    'rs1799884': {'gene': 'GCK', 'pdb': '1V4S', 'drug': 'GCK_Activator'},
    'rs11212617': {'gene': 'ATM', 'pdb': '7SAY', 'drug': 'Metformine'},
    'rs5219': {'gene': 'KCNJ11', 'pdb': '6C3O', 'drug': 'Gliclazide'},
    'rs4124874': {'gene': 'SLC22A1', 'pdb': '8H66', 'drug': 'Bypass_Molecule'},
    'rs7903146': {'gene': 'TCF7L2', 'pdb': '2G4B', 'drug': 'GLP1_Analogue'}
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

def lancer_docking_v2(rsid, info):
    # 1. Préparation du ligand
    ligand = preparer_ligand(info['drug'], MOLECULES_SMILES[info['drug']])
    
    # 2. Vérification des fichiers cibles
    receptor_pdb = f"{info['pdb']}.pdb"
    receptor_pdbqt = f"{info['pdb']}.pdbqt"
    if not os.path.exists(receptor_pdbqt): return "PDBQT manquant"
    
    # 3. INNOVATION : Appel à pocket_finder pour trouver la serrure en temps réel
    print(f"   [AI] Recherche de la serrure profonde pour {info['pdb']}...")
    coords = trouver_poche_profonde(info['pdb'])
    if not coords: return "Erreur PocketFinder"
    cx, cy, cz = coords

    # 4. Nettoyage du récepteur
    nettoyage_vina_pro(receptor_pdbqt)

    # 5. Configuration Vina
    conf = f"conf_{rsid}.txt"
    with open(conf, "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\nligand = {ligand}\n")
        f.write(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n")
        f.write("size_x = 22.0\nsize_y = 22.0\nsize_z = 22.0\nexhaustiveness = 10\n")

    # 6. Exécution
    try:
        print(f"   [VINA] Docking en cours au point {cx:.1f}, {cy:.1f}, {cz:.1f}...")
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        if res.returncode != 0: return "Vina Error"
        
        for line in res.stdout.split('\n'):
            if "   1 " in line:
                return line.split()[1]
    except: return "Erreur Système"
    return "N/A"

def main():
    print("="*90)
    print("  DIABETE-GENE V2 : DOCKING AUTOMATISÉ AVEC RECHERCHE DE SERRURE")
    print("="*90)
    
    if not os.path.exists("patient_expert.vcf"):
        print("Erreur : patient_expert.vcf introuvable.")
        return
        
    reader = vcf.Reader(filename="patient_expert.vcf")
    for record in reader:
        if record.ID in TABLEAU_STRATEGIQUE:
            info = TABLEAU_STRATEGIQUE[record.ID]
            print(f"\nANALYSE MUTATION : {record.ID} ({info['gene']})")
            score = lancer_docking_v2(record.ID, info)
            
            if isinstance(score, str) and "-" in score:
                print(f"   ✅ SUCCÈS ! Affinité : {score} kcal/mol")
            else:
                print(f"   ⚠️ ÉCHEC : {score}")

if __name__ == "__main__":
    main()
