import os, subprocess, re, vcf
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pocketfinderv12 import trouver_poche_profonde 

VINA_EXE = "vina.exe" 
VCF_FILE = "patient_complet.vcf"

# --- BASE DE DONNÉES EXHAUSTIVE (Gène -> PDB -> Médicaments à tester) ---
# On teste plusieurs molécules par gène pour trouver la plus affine
DB_DIABETE = {
    'GCK': {
        'pdb': '1V4S', 
        'drugs': {
            'Activateur_GCK': 'CC1=CC=C(C=C1)S(=O)(=O)N',
            'Dorzagliatin': 'CN1C=C(N=C1)C(=O)N[C@@H](C)C2=CC=C(S(=O)(=O)C)C=C2',
            'Piragliatin': 'CC(C)OC1=C(C=C(C=N1)S(=O)(=O)C2CC2)NC(=O)C3=NN(C=C3C)C'
        }
    },
    'KCNJ11': {
        'pdb': '6C3O',
        'drugs': {
            'Gliclazide': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3',
            'Glibenclamide': 'CCNC1=CC=C(C=C1)S(=O)(=O)NC(=O)NC2CCCCC2',
            'Repaglinide': 'CCCC1=C(C=CC(=C1)CC(=O)O)OCC2=CC=CC=C2N3CCCCC3'
        }
    },
    'INSR': {
        'pdb': '1IRK',
        'drugs': {
            'Insulin_Mimetic': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O',
            'Demethylasterriquinone': 'CC(C)=CCC1=C(C(=C(C(=C1O)CC=C(C)C)C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O)O'
        }
    },
    'SLC22A1': {
        'pdb': '8H66',
        'drugs': {
            'Metformine': 'CN(C)C(=N)N=C(N)N',
            'Berberine': 'COc1c2c(cc3c1OC)C[n+]4c(c2)cc5c(c4CC3)cc6c(c5)OCO6'
        }
    },
    'PPARG': {
        'pdb': '5YCP',
        'drugs': {
            'Pioglitazone': 'CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3',
            'Rosiglitazone': 'CN(C)C1=CC=NC=C1N(C)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3'
        }
    }
}

# --- MAPPING RSID -> GÈNE ---
RSID_MAP = {
    'rs1799884': 'GCK', 'rs11212617': 'ATM', 'rs1801278': 'INSR', 
    'rs5219': 'KCNJ11', 'rs7903146': 'TCF7L2', 'rs4124874': 'SLC22A1',
    'rs1801282': 'PPARG'
}

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

def executer_docking_v7(pdb_id, drug_name, smiles):
    try:
        cx, cy, cz = trouver_poche_profonde(pdb_id)
        ligand = preparer_ligand_v7(drug_name, smiles)
        receptor_raw = f"{pdb_id}.pdbqt"
        if not os.path.exists(receptor_raw): return "PDBQT manquant"
        
        receptor = nettoyer_recepteur_vina(receptor_raw)
        cmd = [VINA_EXE, "--receptor", receptor, "--ligand", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
               "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "8"]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if os.path.exists(receptor): os.remove(receptor)
        scores = re.findall(r"-\d+\.\d+", result.stdout)
        return scores[0] if scores else "N/A"
    except Exception as e: return f"Err: {e}"

def main():
    print("="*80 + "\n  DIABETE-GENE V9 : SCREENING EXHAUSTIF (DRUG REPURPOSING)\n" + "="*80)
    
    # 1. Extraction patient
    vcf_reader = vcf.Reader(filename=VCF_FILE)
    patient_rsids = [record.ID for record in vcf_reader if record.ID]
    
    # 2. Screening
    for rsid in patient_rsids:
        if rsid in RSID_MAP:
            gene = RSID_MAP[rsid]
            if gene in DB_DIABETE:
                target = DB_DIABETE[gene]
                print(f"\n🧬 MUTATION DÉTECTÉE : {rsid} (Gène: {gene})")
                print(f"🔬 Screening des molécules sur la cible {target['pdb']}...")
                
                results = []
                for drug, smiles in target['drugs'].items():
                    score = executer_docking_v7(target['pdb'], drug, smiles)
                    print(f"   [-] {drug.ljust(20)} : {score} kcal/mol")
                    if score != "N/A" and "Err" not in score:
                        results.append((drug, float(score)))
                
                if results:
                    best_drug = min(results, key=lambda x: x[1])
                    print(f"   🏆 MEILLEURE OPTION : {best_drug[0]} ({best_drug[1]} kcal/mol)")

if __name__ == "__main__":
    main()
