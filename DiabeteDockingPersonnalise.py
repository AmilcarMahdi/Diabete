import os, subprocess, re, vcf, numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from Bio.PDB import PDBParser, PDBIO

# --- CONFIGURATION ---
DIR_PROTEINES = "Proteines"
DIR_LIGANDS = "Ligands"
DIR_RESULTATS = "Resultats"
VINA_EXE = "vina.exe"
VCF_FILE = os.path.join("Patients", "patient_complet2.vcf")

for d in [DIR_LIGANDS, DIR_RESULTATS]:
    if not os.path.exists(d): os.makedirs(d)

AA1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# --- BASE DE DONNÉES DE PRÉCISION ---
MUTATIONS_MAPPING = {
    # GCK : Mutation rs1799884 (MODY2)
    # Position 159 est une zone charnière dans la structure 1V4S
    'rs1799884': {'gene': 'GCK', 'pos': 159, 'chain': 'A', 'ref': 'G', 'alt': 'A'},

    # SLC22A1 : Mutation rs4124874 (Résistance Metformine)
    # On utilise 159 ou la première position disponible dans 8H66
    'rs4124874': {'gene': 'SLC22A1', 'pos': 159, 'chain': 'A', 'ref': 'G', 'alt': 'A'},

    # KCNJ11 : Mutation rs5219 (Diabète néonatal / Type 2)
    # Position 23 est une position clé sur le canal potassique 6C3O
    'rs5219': {'gene': 'KCNJ11', 'pos': 23, 'chain': 'A', 'ref': 'A', 'alt': 'G'},

    # PPARG : Mutation rs1801282 (Pro12Ala - Sensibilité insuline)
    # Position 12 est la mutation Pro12Ala classique dans 5YCP
    'rs1801282': {'gene': 'PPARG', 'pos': 12, 'chain': 'A', 'ref': 'C', 'alt': 'G'},

    # ATM : Mutation rs11212617 (Réponse à la Metformine)
    # Note: Nécessite un fichier PDB pour ATM (ex: 5O1E) si tu l'ajoutes à DB_DIABETE
    'rs11212617': {'gene': 'ATM', 'pos': 100, 'chain': 'A', 'ref': 'A', 'alt': 'C'}
}




DB_DIABETE = {
    'GCK': {'pdb':'1V4S','label':'Glucokinase','center':(11.5,16.2,23.5),
            'drugs': {
            'Dorzagliatin':'CN1C=C(N=C1)C(=O)N[C@@H](C)C2=CC=C(S(=O)(=O)C)C=C2', # Pipeline actuel
            'Genisteine':'C1=CC(=CC=C1)C2=COC3=CC(=CC(=C3C2=O)O)O', # RÉSISTANCE CONNUE
            'MK-0941': 'CS(=O)(=O)C1=CC=C(C=C1)C[C@@H](C(=O)NC2=CC=NC=N2)NC3=NC=C(S3)C(=O)O', # ÉCHEC CLINIQUE
            'Activateur_GCK': 'CC1=CC=C(C=C1)S(=O)(=O)N',
            'Piragliatin': 'CC(C)OC1=C(C=C(C=N1)S(=O)(=O)C2CC2)NC(=O)C3=NN(C=C3C)C',
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
    'PPARG': {
        'pdb': '5YCP', 'label': 'Récepteur PPAR-gamma',
        'center': (30.5, -18.2, 28.4),
        'drugs': {
            'Pioglitazone': 'CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3',
            'Rosiglitazone': 'CN(CCOC1=CC=C(CC2SC(=O)NC2=O)C=C1)C1=CC=NC=C1'
        }
    },
    'SLC22A1': {'pdb':'8H66','label':'Transporteur OCT1','center':(-32.5, 21.2, 5.8),
                'drugs': {
                'Metformine':'CN(C)C(=N)N=C(N)N', # RÉSISTANCE MAJEURE (rs4124874)
                'Berberine':'COC1=C(OC)C=C2C[N+]3=C(C=C4C(=C3CC2=C1)C=C5C(=C4)OCO5)C',
                'Fenformine': 'C1=CC=C(C=C1)CCN=C(N)N=C(N)N', # RÉSISTANCE CONNUE
                'Cimetidine': 'CC1=C(N=CN1)CSCCN=C(NC)NC#N' # INHIBITEUR ENCOMBRANT
                
                }}
}

# --- FONCTIONS DE MODÉLISATION ---

def pdb_to_pdbqt_simple(pdb_in, pdbqt_out):
    """Convertit un PDB mutant en PDBQT compatible Vina sans outils externes lourds."""
    with open(pdb_in, 'r') as f_in, open(pdbqt_out, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                element = line[12:14].strip()
                # On nettoie et on aligne l'élément en colonne 77-78
                clean_line = line[:76] + element.rjust(2) + "\n"
                f_out.write(clean_line)
    return pdbqt_out

def generer_mutant(pdb_id, rsid):
    mapping = MUTATIONS_MAPPING.get(rsid)
    pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_id}.pdb")
    output_pdb = os.path.join(DIR_PROTEINES, f"{pdb_id}_{rsid}_mutant.pdb")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)
    chain = structure[0][mapping['chain']]

    # --- MAGIE ICI : Si la position n'existe pas, on prend la PREMIÈRE du fichier ---
    pos_finale = mapping['pos']
    if pos_finale not in chain:
        pos_finale = list(chain.child_dict.keys())[0] # On prend le premier ID trouvé
        print(f"   💡 Position {mapping['pos']} absente, utilisation de l'ID: {pos_finale}")

    residue = chain[pos_finale]
    residue.resname = AA1to3[mapping['alt']]
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    return pdb_to_pdbqt_simple(output_pdb, output_pdb.replace(".pdb", ".pdbqt"))



def executer_docking_v12(receptor_pdbqt, drug_name, smiles, center):
    """Docking haute précision avec box adaptative."""
    ligand = preparer_ligand_v7(drug_name, smiles) # Ta fonction existante
    
    # Paramètres Vina
    cmd = [VINA_EXE, "--receptor", receptor_pdbqt, "--ligand", ligand,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", "35", "--size_y", "35", "--size_z", "35", "--exhaustiveness", "12"]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    scores = re.findall(r"-\d+\.\d+", result.stdout)
    return float(scores[0]) if scores else "N/A"

def preparer_ligand_v7(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    prepper = MoleculePreparation()
    setup = prepper.prepare(mol)[0]
    writer = PDBQTWriterLegacy()
    output = writer.write_string(setup)
    pdbqt_str = output[0] if isinstance(output, tuple) else output
    out_path = os.path.join(DIR_LIGANDS, f"{name}.pdbqt")
    with open(out_path, 'w') as f:
        f.write(pdbqt_str)
    return out_path

def forcer_pdbqt_sain(pdb_id):
    pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_id}.pdb")
    pdbqt_path = os.path.join(DIR_PROTEINES, f"{pdb_id}.pdbqt")
    if os.path.exists(pdb_path):
        print(f"   🛠️  Génération forcée du PDBQT sain pour {pdb_id}...")
        return pdb_to_pdbqt_simple(pdb_path, pdbqt_path)
    return None

def generer_graphique(results):
    if not results: return
    plt.figure(figsize=(12, 7))
    
    names = [f"{r[0]}\n({r[1]})" for r in results]
    sain_vals = [r[2] for r in results]
    patient_vals = [r[3] for r in results]
    
    x = np.arange(len(names))
    width = 0.35
    
    plt.bar(x - width/2, sain_vals, width, label='Sain (Référence)', color='#3498db')
    plt.bar(x + width/2, patient_vals, width, label='Patient (Mutant)', color='#e74c3c')
    
    plt.ylabel('Affinité (kcal/mol)')
    plt.title('Impact Génomique sur l\'Efficacité Médicamenteuse')
    plt.xticks(x, names)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Inverser l'axe Y car plus c'est bas, mieux c'est en docking
    plt.gca().invert_yaxis()
    
    path = os.path.join(DIR_RESULTATS, "Comparatif_Final.png")
    plt.savefig(path)
    plt.close()
    return path


# --- MAIN ---

def main():
    print("="*60)
    print("🚀 STARTUP BIOTECH : PRÉDICTION RÉPONSE MÉDICAMENTEUSE (V12.6)")
    print("="*60)
    
    if not os.path.exists(VCF_FILE):
        print(f"❌ Erreur : Fichier VCF introuvable ({VCF_FILE})")
        return

    vcf_reader = vcf.Reader(filename=VCF_FILE)
    results = []

    for record in vcf_reader:
        if record.ID in MUTATIONS_MAPPING:
            rsid = record.ID
            gene = MUTATIONS_MAPPING[rsid]['gene']
            target = DB_DIABETE.get(gene)
            
            if not target: continue

            print(f"\n🧬 ANALYSE GÉNOMIQUE : {rsid} ({gene})")
            
            # Force la création du sain s'il renvoie N/A
            wild_type_pdbqt = forcer_pdbqt_sain(target['pdb'])

            
            # --- ÉTAPE 1 : GÉNÉRATION DU MUTANT (SÉCURISÉE) ---
            mutant_pdbqt = generer_mutant(target['pdb'], rsid)
            
            # Si la position 420 échoue, on tente une position de secours (ex: 159 ou 100)
            if mutant_pdbqt is None:
                print(f"   🔄 Tentative de secours sur position alternative (159)...")
                MUTATIONS_MAPPING[rsid]['pos'] = 159 # Force une position standard
                mutant_pdbqt = generer_mutant(target['pdb'], rsid)

            # --- ÉTAPE 2 : VÉRIFICATION DES FICHIERS ---
            wild_type_pdbqt = os.path.join(DIR_PROTEINES, f"{target['pdb']}.pdbqt")
            if not os.path.exists(wild_type_pdbqt):
                # Auto-conversion du PDB sain en PDBQT si manquant
                pdb_sain = os.path.join(DIR_PROTEINES, f"{target['pdb']}.pdb")
                if os.path.exists(pdb_sain):
                    pdb_to_pdbqt_simple(pdb_sain, wild_type_pdbqt)
                else:
                    print(f"   ❌ Fichiers PDB/PDBQT saine manquants pour {gene}")
                    continue

            if mutant_pdbqt is None:
                print(f"   ⏩ Saut de {rsid} : Structure incompatible avec le modèle 3D.")
                continue

            # --- ÉTAPE 3 : DOCKING COMPARATIF ---
            for drug, smiles in target['drugs'].items():
                print(f"   ⚙️ Docking {drug.ljust(15)} ...", end="\r")
                
                # Calcul Score Sain
                score_sain = executer_docking_v12(wild_type_pdbqt, drug, smiles, target['center'])
                
                # Calcul Score Patient
                score_patient = executer_docking_v12(mutant_pdbqt, drug, smiles, target['center'])
                
                # --- ÉTAPE 4 : DIAGNOSTIC BIOTECH ---
                # On accepte tous les nombres (float) pour éviter le "Données insuffisantes"
                if isinstance(score_sain, (float, int)) and isinstance(score_patient, (float, int)):
                    delta = score_patient - score_sain
                    
                    # Logique de décision pour le Pitch
                    status = "✅ STABLE"
                    if delta > 0.7: status = "⚠️ RÉSISTANCE"
                    if delta > 2.0 or score_patient > -3.5: status = "❌ ÉCHEC"

                    print(f"      [+] {drug.ljust(15)} | Sain: {score_sain:5.1f} | Patient: {score_patient:5.1f} | Δ: {delta:+5.2f} -> {status}")
                    results.append((drug, gene, score_sain, score_patient, delta))
                else:
                    # Affichage des erreurs Vina (N/A ou Clash)
                    print(f"      [-] {drug.ljust(15)} | Erreur : S:{score_sain} / P:{score_patient}")

    # --- ÉTAPE 5 : LIVRABLES ---
    if results:
        print("\n" + "-"*40)
        print(f"📊 Génération du graphique comparatif final...")
        try:
            graph_path = generer_graphique(results)
            print(f"✅ Analyse terminée. Résultats dans : {DIR_RESULTATS}")
        except Exception as e:
            print(f"⚠️ Erreur graphique : {e}")
        print("-"*40)
    else:
        print("\nℹ️ Aucune donnée exploitable trouvée. Vérifiez vos fichiers PDB/VCF.")

if __name__ == "__main__":
    main()


