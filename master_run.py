import docking_pipeline as dp
import ml_engine as ml
import visualizer as vis
import vcf
import os

# 1. Définition des données
#MUTATIONS = {'rs1799884': {'gene': 'GCK', 'pos': 159, 'chain': 'A', 'ref': 'G', 'alt': 'A'}}
MUTATIONS = {
    # GCK : Mutation rs1799884 (MODY2)
    # Position 159 est une zone charnière dans la structure 1V4S
    'rs1799884': {'gene': 'GCK', 'pos': 159, 'chain': 'A', 'ref': 'G', 'alt': 'A'},

    # SLC22A1 : Mutation rs4124874 (Résistance Metformine)
    # On utilise 159 ou la première position disponible dans 8H66
    'rs4124874': {'gene': 'SLC22A1', 'pos': 500, 'chain': 'A', 'ref': 'G', 'alt': 'A'},

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
#DB = {'GCK': {'pdb':'1V4S','center':(11.5,16.2,23.5), 'drugs': {'Dorzagliatin':'SMILES_ICI'}}}
DB = {
    'GCK': {'pdb':'1V4S','label':'Glucokinase',#'center':(11.5,16.2,23.5),
            'drugs': {
            'Dorzagliatin':'CN1C=C(N=C1)C(=O)N[C@@H](C)C2=CC=C(S(=O)(=O)C)C=C2', # Pipeline actuel
            'Genisteine':'C1=CC(=CC=C1)C2=COC3=CC(=CC(=C3C2=O)O)O', # RÉSISTANCE CONNUE
            'MK-0941': 'CS(=O)(=O)C1=CC=C(C=C1)C[C@@H](C(=O)NC2=CC=NC=N2)NC3=NC=C(S3)C(=O)O' # ÉCHEC CLINIQUE
            """
            'Activateur_GCK': 'CC1=CC=C(C=C1)S(=O)(=O)N',
            'Piragliatin': 'CC(C)OC1=C(C=C(C=N1)S(=O)(=O)C2CC2)NC(=O)C3=NN(C=C3C)C',
            'Quercetine': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C2=O)C=CC(=C3)O)O)O)O',
            'Resveratrol': 'C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O',
            'Apigenine': 'C1=CC(=CC=C1)C2=CC(=O)C3=C(O2)C=C(C=C3O)O',
            'Kaempferol': 'C1=CC(=CC=C1)C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O',
            'Myricetine': 'C1=CC(=C(C(=C1)O)O)C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O'
            """
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
    'SLC22A1': {'pdb':'8H66','label':'Transporteur OCT1',
                'drugs': {
                'Metformine':'CN(C)C(=N)N=C(N)N', # RÉSISTANCE MAJEURE (rs4124874)
                'Berberine':'COC1=C(OC)C=C2C[N+]3=C(C=C4C(=C3CC2=C1)C=C5C(=C4)OCO5)C',
                'Fenformine': 'C1=CC=C(C=C1)CCN=C(N)N=C(N)N', # RÉSISTANCE CONNUE
                'Cimetidine': 'CC1=C(N=CN1)CSCCN=C(NC)NC#N' # INHIBITEUR ENCOMBRANT
                
        }
    },
    'ATM': {'pdb': '5O1E',
            'drugs': {'Metformine': 'CN(C)C(=N)N=C(N)N', 'KU-55933': 'Cc1cc(O)c2c(c1)oc(-c1cccc(CN3CCOCC3)c1)cc2=O'}}
}


def extraire_mutations_patient(vcf_path, mapping):
    """Lit le VCF et retourne les mutations correspondantes au mapping."""
    detectees = {}
    if not os.path.exists(vcf_path):
        print(f"⚠️ Fichier VCF introuvable : {vcf_path}. Utilisation de tout le mapping pour test.")
        return mapping # Retourne tout pour éviter de bloquer si le fichier manque
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split('\t')
            rsid = parts[2] # Colonne ID du VCF
            if rsid in mapping:
                print(f"   ✅ Mutation trouvée chez le patient : {rsid}")
                detectees[rsid] = mapping[rsid]
    return detectees

# --- 3. EXECUTION PRINCIPALE ---

if __name__ == "__main__":
    print("\n" + "="*40)
    print("🚀 STARTUP BIOTECH : PIPELINE V14.2")
    print("="*40)

    # Chemin vers ton VCF
    vcf_patient = os.path.join("Patients", "patient_complet2.vcf")

    # ÉTAPE 0 : Analyse Génomique
    mutations_patient = extraire_mutations_patient(vcf_patient, MUTATIONS)

    if not mutations_patient:
        print("❌ Aucune mutation d'intérêt détectée.")
    else:
        # ÉTAPE 1 : Docking Moléculaire (Physique)
        print("\n--- ÉTAPE 1 : DOCKING MULTIPLE ---")
        dp.generate_raw_results(mutations_patient, DB)

        # ÉTAPE 2 : Machine Learning (Feature Engineering)
        print("\n--- ÉTAPE 2 : MACHINE LEARNING ---")
        ml.train_model()

        # ÉTAPE 3 : Visualisation (Rapport Q1)
        print("\n--- ÉTAPE 3 : GÉNÉRATION DES FIGURES ---")
        vis.produce_figures()

    print("\n" + "="*40)
    print("✨ PIPELINE TERMINÉ AVEC SUCCÈS")
    print("="*40)
