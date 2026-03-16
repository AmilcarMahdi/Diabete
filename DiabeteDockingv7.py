import os, subprocess, re, vcf
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pocketfinderv12 import trouver_poche_profonde 

VINA_EXE = "vina.exe" 
VCF_FILE = "patient_complet.vcf" # <-- Mets le nom de ton fichier VCF ici

STRATEGIE = {
    'rs1801278': {
        'gene': 'INSR', 'pdb': '1IRK', 'drug': 'Insulin_Mimetic', 
        'smiles': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O'
    },
    'rs5219': {
        'gene': 'KCNJ11', 'pdb': '6C3O', 'drug': 'Gliclazide', 
        'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3'
    },
    'rs1799884': {
        'gene': 'GCK', 'pdb': '1V4S', 'drug': 'Activateur_GCK', 
        'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N'
    }
    # Ajoute les autres ici...
}

def extraire_mutations_patient(vcf_path):
    """Lit le VCF et retourne la liste des RSIDs trouvés."""
    rsids_patient = []
    if not os.path.exists(vcf_path):
        print(f"⚠️ Fichier {vcf_path} introuvable.")
        return []
    
    vcf_reader = vcf.Reader(filename=vcf_path)
    for record in vcf_reader:
        if record.ID:
            rsids_patient.append(record.ID)
    return rsids_patient

def nettoyer_recepteur_vina(pdbqt_file):
    temp_file = pdbqt_file.replace(".pdbqt", "_clean.pdbqt")
    types_valides = {'C', 'A', 'N', 'NA', 'OA', 'SA', 'S', 'P', 'H', 'HD', 'F', 'CL', 'BR', 'I'}
    with open(pdbqt_file, 'r') as f_in, open(temp_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                line = line.replace("N +", "N  ").replace("O -", "O  ").replace("C +", "C  ")
                element = line[76:78].strip().upper()
                if element in types_valides: f_out.write(line)
            elif line.startswith(("TER", "ENDMDL")): f_out.write(line)
    return temp_file

def preparer_ligand_v7(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt_str, _, _ = PDBQTWriterLegacy().write_string(setup)
    out = f"{name}.pdbqt"
    with open(out, "w") as f:
        for line in pdbqt_str.split('\n'):
            if line.startswith(("ATOM", "HETATM")):
                line = line[:77] + line[12:14].strip().ljust(2) + line[79:]
            f.write(line + "\n")
    return out

def executer_docking_v7(rsid, info):
    try:
        cx, cy, cz = trouver_poche_profonde(info['pdb'])
        ligand = preparer_ligand_v7(info['drug'], info['smiles'])
        receptor_raw = f"{info['pdb']}.pdbqt"
        if not os.path.exists(receptor_raw): return "PDBQT manquant"
        
        receptor = nettoyer_recepteur_vina(receptor_raw)
        cmd = [VINA_EXE, "--receptor", receptor, "--ligand", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
               "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "8"]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if os.path.exists(receptor): os.remove(receptor)
        scores = re.findall(r"-\d+\.\d+", result.stdout)
        return scores if scores else "N/A"
    except Exception as e: return f"Erreur: {e}"

def main():
    print("="*80 + "\n  DIABETE-GENE V8 : ANALYSE BASÉE SUR LE VCF PATIENT\n" + "="*80)
    
    # 1. Scanner le patient
    mutations_trouvees = extraire_mutations_patient(VCF_FILE)
    print(f"🧬 Mutations détectées dans le VCF : {len(mutations_trouvees)}")
    
    # 2. Filtrer et calculer
    matches = [m for m in mutations_trouvees if m in STRATEGIE]
    
    if not matches:
        print("✅ Aucune mutation à risque (parmi la liste cible) n'a été trouvée.")
        return

    print(f"🎯 {len(matches)} mutation(s) correspondent à notre base de données. Lancement du docking...\n")
    
    for rsid in matches:
        info = STRATEGIE[rsid]
        print(f"[CIBLE TROUVÉE] {rsid} (Gène: {info['gene']})")
        score = executer_docking_v7(rsid, info)
        print(f"   >>> Médicament : {info['drug']} | Affinité : {score} kcal/mol")

if __name__ == "__main__":
    main()
