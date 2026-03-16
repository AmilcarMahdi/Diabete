import os, subprocess, re, vcf
import matplotlib.pyplot as plt
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pocketfinderv12 import trouver_poche_profonde 

# --- CONFIGURATION ---
VINA_EXE = "vina.exe" 
VCF_FILE = "./Patients/patient_complet.vcf"

# --- BASE DE DONNÉES CIBLES & DROGUES ---
DB_DIABETE = {
    'GCK': {
        'pdb': '1V4S', 'label': 'Glucokinase',
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
        'drugs': {
            'Gliclazide': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3',
            'Glibenclamide': 'ClC1=CC(=C(C=C1)C(=O)NCC2=CC=C(S(=O)(=O)NC(=O)NC3CCCCC3)C=C2)OC',
            'Repaglinide': 'CCCC1=C(C=CC(=C1)CC(=O)O)OCC2=CC=CC=C2N3CCCCC3',
            'Verapamil': 'CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC'
        }
    },
    'SLC22A1': {
        'pdb': '8H66', 'label': 'Transporteur OCT1',
        'drugs': {
            'Metformine': 'CN(C)C(=N)N=C(N)N',
            'Berberine': 'COC1=C(OC)C=C2C[N+]3=C(C=C4C(=C3CC2=C1)C=C5C(=C4)OCO5)C'
        }
    },
    'PPARG': {
        'pdb': '5YCP', 'label': 'Récepteur PPAR-gamma',
        'drugs': {
            'Pioglitazone': 'CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3',
            'Rosiglitazone': 'CN(CCOC1=CC=C(CC2SC(=O)NC2=O)C=C1)C1=CC=NC=C1'
        }
    }
}

# --- MAPPING MUTATIONS ---
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
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt_str = PDBQTWriterLegacy().write_string(setup)[0]
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
        if not os.path.exists(receptor_raw): return "Missing PDB"
        
        receptor = nettoyer_recepteur_vina(receptor_raw)
        # Augmentation de la box à 30 pour éviter les collisions Metformine
        cmd = [VINA_EXE, "--receptor", receptor, "--ligand", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
               "--size_x", "30", "--size_y", "30", "--size_z", "30", "--exhaustiveness", "8"]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if os.path.exists(receptor): os.remove(receptor)
        scores = re.findall(r"-\d+\.\d+", result.stdout)
        return float(scores[0]) if scores else "N/A"
    except Exception as e: return f"Err: {e}"

def generer_graphique(results_final):
    plt.figure(figsize=(10, 6))
    drugs = [f"{r[0]}\n({r[1]})" for r in results_final]
    scores = [r[2] for r in results_final]
    colors = ['green' if s < -6 else 'orange' for s in scores]
    
    plt.bar(drugs, scores, color=colors)
    plt.axhline(y=-6.0, color='red', linestyle='--', label='Seuil Efficacité (-6)')
    plt.ylabel("Affinité (kcal/mol)")
    plt.title("Comparaison des Affinités Médicamenteuses par Mutation")
    plt.gca().invert_yaxis() # Inverser car plus c'est négatif, mieux c'est
    plt.savefig("resultats_docking.png")
    plt.close()

def main():
    print("="*60 + "\n  V10 : GENOTYPAGE ET REPURPOSING PHARMACOLOGIQUE\n" + "="*60)
    
    if not os.path.exists(VCF_FILE): return print("Fichier VCF absent.")
    vcf_reader = vcf.Reader(filename=VCF_FILE)
    patient_mutations = [r.ID for r in vcf_reader if r.ID in MUTATIONS_CIBLES]
    
    rapport_data = []
    
    for rsid in patient_mutations:
        gene = MUTATIONS_CIBLES[rsid]
        if gene in DB_DIABETE:
            target = DB_DIABETE[gene]
            print(f"\n🧬 Mutation : {rsid} | Cible : {target['label']}")
            
            for drug, smiles in target['drugs'].items():
                score = executer_docking_v7(target['pdb'], drug, smiles)
                
                # Filtrage scores aberrants
                if isinstance(score, float) and score > -15.0:
                    print(f"   [-] {drug.ljust(15)} : {score} kcal/mol")
                    rapport_data.append((drug, gene, score))
                else:
                    print(f"   [-] {drug.ljust(15)} : Ignoré (Score: {score})")

    if rapport_data:
        # 1. Graphique
        generer_graphique(rapport_data)
        
        # 2. PDF
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("helvetica", 'B', 16)
        pdf.cell(200, 10, "Rapport d'Analyse Pharmaco-Genomique", 
                 new_x="LMARGIN", new_y="NEXT", align='C')
        pdf.set_font("helvetica", size=12)
        pdf.ln(10)
        pdf.cell(200, 10, f"Patient VCF : {VCF_FILE}", 
                 new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        
        for drug, gene, score in rapport_data:
            statut = "EFFICACE" if score < -6 else "FAIBLE"
            pdf.cell(200, 10, f"Gene: {gene} | Drug: {drug} | Score: {score} kcal/mol | [{statut}]", 
                     new_x="LMARGIN", new_y="NEXT")
        
        pdf.image("resultats_docking.png", x=10, y=100, w=180)
        pdf.output("Rapport_Diabete_V10.pdf")
        print("\n✅ Analyse terminée. Rapport PDF et Graphique générés.")

if __name__ == "__main__":
    main()
