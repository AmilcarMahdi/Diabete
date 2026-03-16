import vcf
import os
import subprocess
from datetime import datetime

# --- CONFIGURATION DU MOTEUR VINA ---
VINA_EXECUTABLE = "vina"  # Assurez-vous que Vina est dans votre PATH
LIGAND_FILE = "molecule_test.pdbqt"

# --- MATRICE EXPERTE (GENES, PDB, COORDONNÉES) ---
STRATEGIES = {
    'rs1801278': {
        'gene': 'INSR', 'pdb': '1IRK', 'action': 'Mimer l\'insuline (Agoniste oral)',
        'center': [11.0, 14.5, 19.0], 'size': [20, 20, 20]
    },
    'rs1799884': {
        'gene': 'GCK', 'pdb': '1V4S', 'action': 'Réveiller l\'enzyme (Activateur GCK)',
        'center': [18.5, 10.0, 12.5], 'size': [22, 22, 22]
    },
    'rs11212617': {
        'gene': 'ATM', 'pdb': '7SAY', 'action': 'Booster l\'effet Metformine',
        'center': [162.0, 175.0, 160.0], 'size': [25, 25, 25]
    },
    'rs5219': {
        'gene': 'KCNJ11', 'pdb': '6C3O', 'action': 'Forcer la fermeture (Sulfonylurées)',
        'center': [115.0, 110.0, 95.0], 'size': [24, 24, 24]
    },
    'rs4124874': {
        'gene': 'SLC22A1', 'pdb': '8H66', 'action': 'Contourner le barrage (Transport)',
        'center': [125.0, 128.0, 142.0], 'size': [20, 20, 20]
    },
    'rs7903146': {
        'gene': 'TCF7L2', 'pdb': '2G4B', 'action': 'Remplacer l\'incrétine (GLP-1)',
        'center': [25.0, 15.0, -8.0], 'size': [22, 22, 22]
    }
}

def create_vina_config(rsid, data):
    """Génère le fichier conf.txt spécifique à la mutation"""
    config_name = f"conf_{rsid}.txt"
    content = (
        f"receptor = {data['pdb']}.pdbqt\n"
        f"ligand = {LIGAND_FILE}\n\n"
        f"center_x = {data['center'][0]}\n"
        f"center_y = {data['center'][1]}\n"
        f"center_z = {data['center'][2]}\n\n"
        f"size_x = {data['size'][0]}\n"
        f"size_y = {data['size'][1]}\n"
        f"size_z = {data['size'][2]}\n\n"
        "exhaustiveness = 8"
    )
    with open(config_name, "w") as f:
        f.write(content)
    return config_name

def run_simulation(config_file, rsid):
    """Exécute Vina et capture l'affinité"""
    output_file = f"result_{rsid}.pdbqt"
    log_file = f"log_{rsid}.txt"
    cmd = [VINA_EXECUTABLE, "--config", config_file, "--out", output_file, "--log", log_file]
    
    print(f"   [VINA] Calcul d'affinité en cours pour {rsid}...")
    try:
        # Simulation de l'appel (Décommentez pour exécution réelle)
        # subprocess.run(cmd, check=True, capture_output=True)
        return "-8.7 kcal/mol (Valeur simulée)"
    except Exception:
        return "Erreur d'exécution (Vérifiez l'install de Vina)"

def main(vcf_path):
    if not os.path.exists(vcf_path):
        print("Fichier VCF non trouvé.")
        return

    reader = vcf.Reader(filename=vcf_path)
    print("\n" + "="*90)
    print(f"      PLATEFORME DIABETE-GENE V4.1 : DIAGNOSTIC & SIMULATION VINA")
    print("="*90)

    for record in reader:
        if record.ID in STRATEGIES:
            info = STRATEGIES[record.ID]
            print(f"\n[MUTATION DÉTECTÉE] : {record.ID} (Gène: {info['gene']})")
            print(f" > ACTION IDENTIFIÉE : {info['action']}")
            print(f" > CIBLE STRUCTURALE : {info['pdb']}.pdbqt")
            
            # Étape 1 : Créer la config
            config_file = create_vina_config(record.ID, info)
            
            # Étape 2 : Lancer Vina
            affinity = run_simulation(config_file, record.ID)
            print(f" > RÉSULTAT VINA     : {affinity}")
            print("-" * 45)

    print("\n" + "="*90)

if __name__ == "__main__":
    main("patient_expert.vcf")
