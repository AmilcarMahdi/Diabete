import os, subprocess, re, vcf
import matplotlib.pyplot as plt
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pocketfinderv13 import trouver_poche_profonde 

# --- CONFIGURATION DES CHEMINS ---
DIR_PATIENTS = "Patients"
DIR_PROTEINES = "Proteines"
DIR_LIGANDS = "Ligands"
DIR_RESULTATS = "Resultats"

# Création des dossiers si absents
for d in [DIR_LIGANDS, DIR_RESULTATS]:
    if not os.path.exists(d): os.makedirs(d)

VINA_EXE = "vina.exe" 
VCF_FILE = os.path.join(DIR_PATIENTS, "patient_complet.vcf")

# --- BASE DE DONNÉES CIBLES & DROGUES ---
DB_DIABETE = {
    'GCK': {
        'pdb': '1V4S', 'label': 'Glucokinase',
        'center': (11.5, 16.2, 23.5),
        'drugs': {
            'Activateur_GCK': 'CC1=CC=C(C=C1)S(=O)(=O)N',
            'Dorzagliatin': 'CN1C=C(N=C1)C(=O)N[C@@H](C)C2=CC=C(S(=O)(=O)C)C=C2',
            'Piragliatin': 'CC(C)OC1=C(C=C(C=N1)S(=O)(=O)C2CC2)NC(=O)C3=NN(C=C3C)C',
            'Genisteine': 'C1=CC(=CC=C1)C2=COC3=CC(=CC(=C3C2=O)O)O',
            'Quercetine': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C2=O)C=CC(=C3)O)O)O)O',
            'Resveratrol': 'C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O',
            'Apigenine': 'C1=CC(=CC=C1)C2=CC(=O)C3=C(O2)C=C(C=C3O)O',
            'Kaempferol': 'C1=CC(=CC=C1)C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O',
            'Myricetine': 'C1=CC(=C(C(=C1)O)O)C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O'
        }
    },
    'KCNJ11': {
        'pdb': '6C3O', 'label': 'Canal SUR1',
        'center': (148.2, 128.5, 238.1),
        'drugs': {
            'Gliclazide': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3',
            'Glibenclamide': 'ClC1=CC(=C(C=C1)C(=O)NCC2=CC=C(S(=O)(=O)NC(=O)NC3CCCCC3)C=C2)OC',
            'Repaglinide': 'CCCC1=C(C=CC(=C1)CC(=O)O)OCC2=CC=CC=C2N3CCCCC3',
            'Verapamil': 'CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC'
        }
    },
    'SLC22A1': {
        'pdb': '8H66', 'label': 'Transporteur OCT1',
        'center': (125.4, 115.8, 140.2),
        'drugs': {
            'Metformine': 'CN(C)C(=N)N=C(N)N',
            'Berberine': 'COC1=C(OC)C=C2C[N+]3=C(C=C4C(=C3CC2=C1)C=C5C(=C4)OCO5)C'
        }
    },
    'PPARG': {
        'pdb': '5YCP', 'label': 'Récepteur PPAR-gamma',
        'center': (30.5, -18.2, 28.4),
        'drugs': {
            'Pioglitazone': 'CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3',
            'Rosiglitazone': 'CN(CCOC1=CC=C(CC2SC(=O)NC2=O)C=C1)C1=CC=NC=C1'
        }
    }
}

MUTATIONS_CIBLES = {
    'rs1801278': 'INSR', 'rs1799884': 'GCK', 'rs11212617': 'ATM',
    'rs5219': 'KCNJ11', 'rs4124874': 'SLC22A1', 'rs7903146': 'TCF7L2',
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        prepper = MoleculePreparation()
        setup = prepper.prepare(mol)[0] # Extraction de la LISTE
        
        writer = PDBQTWriterLegacy()
        # --- CORRECTIF CRUCIAL ICI ---
        output = writer.write_string(setup)
        pdbqt_str = output[0] if isinstance(output, tuple) else output
        
        out_path = os.path.join(DIR_LIGANDS, f"{name}.pdbqt")
        
        with open(out_path, "w") as f:
            for line in pdbqt_str.split('\n'):
                if line.startswith(("ATOM", "HETATM")):
                    # Forçage du type d'atome en fin de ligne (colonnes 77-78)
                    atom_type = line[12:14].strip()
                    line = line[:76] + atom_type.rjust(2)
                f.write(line + "\n")
        
        # Debug : Vérifier si le fichier n'est pas vide
        if os.path.getsize(out_path) < 100:
            print(f"   ⚠️ ATTENTION : Fichier {name}.pdbqt vide !")
            
        return out_path
    except Exception as e:
        print(f"   ❌ Erreur préparation {name}: {e}")
        return None



def executer_docking_v7(target_dict, drug_name, smiles):
    try:
        pdb_id = target_dict['pdb']
        receptor_raw = os.path.join(DIR_PROTEINES, f"{pdb_id}.pdbqt")
        
        # 1. Calcul du centre de masse initial
        coords = []
        with open(receptor_raw, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        
        import numpy as np
        atoms = np.array(coords)
        cx, cy, cz = np.mean(atoms, axis=0)

        # 2. SÉCURITÉ SPÉCIFIQUE POUR SLC22A1 (OCT1)
        # Cette protéine est un tunnel ; le centre de masse est souvent dans une paroi.
        if pdb_id == "8H66":
            cx += 5.0  # On décale de 5 Angströms pour sortir du "mur"
            cz += 5.0

        ligand = preparer_ligand_v7(drug_name, smiles)
        receptor_clean = nettoyer_recepteur_vina(receptor_raw)
        
        # 3. Lancement Vina avec Box large (40)
        cmd = [VINA_EXE, "--receptor", receptor_clean, "--ligand", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
               "--size_x", "40", "--size_y", "40", "--size_z", "40", 
               "--exhaustiveness", "20"]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if os.path.exists(receptor_clean): os.remove(receptor_clean)
        
        scores = re.findall(r"-\d+\.\d+", result.stdout)
        final_score = float(scores[0]) if scores else "N/A"

        # 4. FILTRE DE RÉSULTAT : Si c'est un clash, on le marque comme N/A
        if isinstance(final_score, float) and final_score < -15.0:
            return "Clash (Poche pleine)"
            
        return final_score

    except Exception as e: 
        return f"Err: {e}"






def generer_graphique(results_final):
    plt.figure(figsize=(10, 6))
    drugs = [f"{r[0]}\n({r[1]})" for r in results_final]
    scores = [r[2] for r in results_final]
    colors = ['green' if s < -6 else 'orange' for s in scores]
    
    plt.bar(drugs, scores, color=colors)
    plt.axhline(y=-6.0, color='red', linestyle='--', label='Seuil Efficacité (-6)')
    plt.ylabel("Affinité (kcal/mol)")
    plt.title("Comparaison des Affinités Médicamenteuses par Mutation")
    plt.gca().invert_yaxis()
    
    graph_path = os.path.join(DIR_RESULTATS, "resultats_docking.png")
    plt.savefig(graph_path)
    plt.close()
    return graph_path

def main():
    print("="*60 + "\n  ANALYSE DIABÈTE : COORDONNÉES FIXES\n" + "="*60)
    
    if not os.path.exists(VCF_FILE): return print("VCF absent.")
    vcf_reader = vcf.Reader(filename=VCF_FILE)
    
    rapport_data = []
    
    for record in vcf_reader:
        if record.ID in MUTATIONS_CIBLES:
            rsid = record.ID
            gene_name = MUTATIONS_CIBLES[rsid]
            
            if gene_name in DB_DIABETE:
                target = DB_DIABETE[gene_name] # On récupère le dictionnaire du gène
                print(f"\n🧬 Mutation : {rsid} | Gène : {gene_name}")
                
                for drug, smiles in target['drugs'].items():
                    # --- ICI : On passe 'target' (le dictionnaire) et non 'gene_name' ---
                    score = executer_docking_v7(target, drug, smiles)
                    
                    if isinstance(score, float):
                        print(f"   [+] {drug.ljust(15)} : {score} kcal/mol")
                        rapport_data.append((drug, gene_name, score))
                    else:
                        print(f"   [-] {drug.ljust(15)} : {score}")

    # (Suite du code pour graphique et PDF...)


if __name__ == "__main__":
    main()
