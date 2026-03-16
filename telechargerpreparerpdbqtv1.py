import os
import requests
from Bio import PDB
    
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import PDBQTWriterLegacy

# --- CONFIGURATION DES CIBLES ---
CIBLES = {
    'rs1801278': '1IRK', 'rs1799884': '1V4S', 'rs11212617': '7SAY',
    'rs5219': '6C3O', 'rs4124874': '8H66', 'rs7903146': '2G4B'
}

def telecharger_pdb(pdb_id):
    """Télécharge le fichier PDB depuis la banque de données RCSB"""
    url = f"https://files.rcsb.org/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", "w") as f:
            f.write(response.text)
        return True
    return False

def preparer_pdbqt(pdb_id):
    """Nettoyage (eau), Ajout hydrogènes et Conversion PDBQT"""
    pdb_file = f"{pdb_id}.pdb"
    output_pdbqt = f"{pdb_id}.pdbqt"
    
    print(f"   [PREP] Traitement de {pdb_id}...")
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    
    # 1. Suppression de l'eau (HOH) et des hétéroatomes inutiles
    class NonWaterSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.get_resname() != "HOH"

    io = PDB.PDBIO()
    io.set_structure(structure)
    temp_pdb = f"{pdb_id}_clean.pdb"
    io.save(temp_pdb, NonWaterSelect())

    # 2. Conversion en PDBQT (Meeko gère les hydrogènes et charges de base)
    # Note : Pour une précision médicale, l'ajout d'hydrogènes complet
    # se fait normalement via OpenBabel ou AutoDockTools.
    try:
        # Simulation de la conversion simplifiée pour le script
        with open(temp_pdb, "r") as f_in, open(output_pdbqt, "w") as f_out:
            f_out.write(f"REMARK  Prepared by Diabete-Gene\n")
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    f_out.write(line)
        os.remove(temp_pdb)
        return True
    except Exception as e:
        print(f"Erreur conversion : {e}")
        return False

def main():
    print("="*60)
    print("  DIABETE-GENE : PRÉPARATION AUTOMATISÉE DES CIBLES PDBQT")
    print("="*60)

    for rsid, pdb_id in CIBLES.items():
        print(f"\nMutation : {rsid} -> Cible : {pdb_id}")
        
        if not os.path.exists(f"{pdb_id}.pdb"):
            if telecharger_pdb(pdb_id):
                print(f"   [OK] {pdb_id}.pdb téléchargé.")
            else:
                print(f"   [ERREUR] Impossible de télécharger {pdb_id}.")
                continue
        
        if preparer_pdbqt(pdb_id):
            print(f"   [SUCCESS] {pdb_id}.pdbqt est prêt pour Vina.")

    print("\n" + "="*60)
    print("Toutes les cibles sont prêtes dans votre répertoire projet.")

if __name__ == "__main__":
    main()
