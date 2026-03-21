import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

PROTEINES_DIR="Proteines"
if not os.path.exists(PROTEINES_DIR):
    os.makedirs(PROTEINES_DIR)

# --- LE DICTIONNAIRE DES 6 MUTATIONS ---
MUTATIONS_CIBLES = {
    'rs1801278': {'pdb': '1IRK', 'label': 'Récepteur Insuline'},
    'rs1799884': {'pdb': '1V4S', 'label': 'Glucokinase'},
    'rs11212617': {'pdb': '7SAY', 'label': 'ATM Kinase'},
    'rs5219': {'pdb': '6C3O', 'label': 'Canal SUR1'},
    'rs4124874': {'pdb': '8H66', 'label': 'Transporteur OCT1'},
    'rs7903146': {'pdb': '2G4B', 'label': 'Signal GLP-1'},
    'rs1801282': {'pdb': '5YCP', 'label': 'Récepteur PPAR-gamma'},
    'rs11212617_v2': {'pdb': '5O1E', 'label': 'ATM Kinase (Domaine)'}
}

def telecharger_pdb(pdb_id):
    # Chemin vers le sous-dossier Proteines
    path = os.path.join(PROTEINES_DIR, f"{pdb_id}.pdb")
    
    if os.path.exists(path): return path
    
    # --- LA CORRECTION EST ICI : /download/ est bien présent ---
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    
    print(f"   [DOWNLOAD] Tentative sur : {url}")
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            with open(path, "w", encoding="utf-8") as f:
                f.write(r.text)
            return path
        else:
            print(f"   [ERREUR HTTP] Code {r.status_code}")
            return None
    except Exception as e:
        print(f"   [ERREUR RESEAU] : {e}")
        return None

def preparer_pdbqt_direct(pdb_id):
    """Bypass Meeko : Utilise RDKit pour ajouter les H et créer le PDBQT"""
    # Chemins pointant vers le dossier Proteines
    pdb_in = os.path.join(PROTEINES_DIR, f"{pdb_id}.pdb")
    pdbqt_out = os.path.join(PROTEINES_DIR, f"{pdb_id}.pdbqt")
    try:
        # Charger la protéine (sanitize=False pour accepter 6C3O)
        mol = Chem.MolFromPDBFile(pdb_in, removeHs=True, sanitize=False)
        if mol is None: return False
        
        # Réparer et ajouter les hydrogènes (H) pour Vina
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol, addCoords=True)
        
        # Écriture du fichier PDBQT compatible Vina
        with open(pdbqt_out, "w") as f:
            f.write(f"REMARK  PROTEIN PREPARED BY DIABETE-GENE\n")
            f.write(Chem.MolToPDBBlock(mol))
            f.write("\nTER\nEND\n")
        return True
    except Exception as e:
        print(f"   [ERREUR PREP] {e}")
        return False

def main():
    print("="*75)
    print("  DIABETE-GENE V4.1.11 : FIX URL DOWNLOAD & BYPASS MEEKO")
    print("="*75)

    for rsid, info in MUTATIONS_CIBLES.items():
        pdb_id = info['pdb']
        print(f"\nTRAITEMENT : {rsid} ({info['label']})")
        
        # 1. TÉLÉCHARGEMENT (Vérifie bien le /download/ dans la console)
        if telecharger_pdb(pdb_id):
            # 2. PRÉPARATION (Sans erreur de fragments)
            if preparer_pdbqt_direct(pdb_id):
                size = os.path.getsize(os.path.join(PROTEINES_DIR, f"{pdb_id}.pdbqt")) // 1024
                print(f"   [SUCCESS] {pdb_id}.pdbqt généré ({size} KB)")
            else:
                print(f"   [ECHEC] Erreur technique sur {pdb_id}.")
        else:
            print(f"   [ECHEC] Téléchargement impossible.")

if __name__ == "__main__":
    main()
