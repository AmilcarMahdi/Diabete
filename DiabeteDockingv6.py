import vcf, os, subprocess, warnings, re
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pocketfinderv12 import trouver_poche_profonde 

VINA_EXE = "vina.exe" 

STRATEGIE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'panne': 'Binding',    'calcul': True,  'drug': 'Insulin_Mimetic', 'smiles': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O'},
    'rs5219':    {'gene': 'KCNJ11','pdb': '6C3O', 'panne': 'Binding',    'calcul': True,  'drug': 'Gliclazide',      'smiles': 'CC1=CC-C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3'}
}

def nettoyer_recepteur_vina(pdbqt_file):
    """Filtre la protéine et nettoie les charges mal formatées comme 'N +'."""
    temp_file = pdbqt_file.replace(".pdbqt", "_clean.pdbqt")
    
    # Types d'atomes acceptés par AutoDock Vina
    types_valides = {'C', 'A', 'N', 'NA', 'OA', 'SA', 'S', 'P', 'H', 'HD', 'F', 'CL', 'BR', 'I'}
    
    with open(pdbqt_file, 'r') as f_in, open(temp_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                # 1. Nettoyage immédiat des charges collées au type d'atome (ex: 'N +' -> 'N  ')
                line = line.replace("N +", "N  ").replace("O -", "O  ").replace("C +", "C  ")
                
                # 2. Extraction du type d'atome (colonne 77-78)
                element = line[76:78].strip().upper()
                
                # 3. Si l'élément est valide, on écrit
                if element in types_valides:
                    f_out.write(line)
            
            elif line.startswith(("TER", "ENDMDL")):
                f_out.write(line)
                
    return temp_file



def preparer_ligand_v7(name, smiles):
    out = f"{name}.pdbqt"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    prepper = MoleculePreparation()
    setup_list = prepper.prepare(mol)
    setup = setup_list[0]
    writer = PDBQTWriterLegacy()
    pdbqt_str, _, _ = writer.write_string(setup)
    
    with open(out, "w") as f:
        for line in pdbqt_str.split('\n'):
            if line.startswith(("ATOM", "HETATM")):
                line = line[:77] + line[12:14].strip().ljust(2) + line[79:]
            f.write(line + "\n")
    return out

def executer_docking_v7(rsid, info):
    center = trouver_poche_profonde(info['pdb'])
    cx, cy, cz = center
    
    ligand = preparer_ligand_v7(info['drug'], info['smiles'])
    receptor_raw = f"{info['pdb']}.pdbqt"
    
    if not os.path.exists(receptor_raw):
        print(f"   ❌ {receptor_raw} manquant.")
        return "N/A"

    # NETTOYAGE DU RÉCEPTEUR ICI
    receptor = nettoyer_recepteur_vina(receptor_raw)
    
    print(f"   🚀 Vina sur {info['pdb']} (Nettoyé)...")
    
    cmd = [
        VINA_EXE, "--receptor", receptor, "--ligand", ligand,
        "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", "8"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Ménage
    if os.path.exists(receptor): os.remove(receptor)
    
    if result.returncode != 0:
        print(f"   ❌ ERREUR : {result.stderr.strip()[:100]}")
        return "N/A"
    
    scores = re.findall(r"-\d+\.\d+", result.stdout)
    return scores[0] if scores else "N/A"

def main():
    print("="*80 + "\n  DIABETE-GENE V7 : FIX PARSING PDBQT\n" + "="*80)
    for rsid, info in STRATEGIE.items():
        if info['calcul']:
            print(f"[MUTATION] {rsid} ({info['gene']})")
            score = executer_docking_v7(rsid, info)
            print(f"   >>> Affinité : {score} kcal/mol")

if __name__ == "__main__":
    main()
