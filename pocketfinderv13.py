import numpy as np
import os

# --- CONFIGURATION DE L'ARBORESCENCE ---
DIR_PROTEINES = "Proteines"

def trouver_poche_profonde(pdb_input):
    if os.path.exists(pdb_input):
        pdb_path = pdb_input
    else:
        pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_input}.pdbqt")
        if not os.path.exists(pdb_path):
            pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_input}.pdb")

    if not os.path.exists(pdb_path):
        print(f"   [PocketFinderV13] Erreur : Fichier introuvable -> {pdb_path}")
        return None
        
    coords = []
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        if not coords:
            print(f"   [PocketFinderV13] Erreur : Aucune coordonnee dans {pdb_path}")
            return None
            
        atoms = np.array(coords)
        np.random.seed(42) 
        
        min_c = np.min(atoms, axis=0)
        max_c = np.max(atoms, axis=0)
        
        grid = np.random.uniform(min_c, max_c, (5000, 3)) 
        
        meilleure_poche = np.mean(atoms, axis=0)
        max_voisins = 0 
        
        for point in grid:
            dist = np.linalg.norm(atoms - point, axis=1)
            
            if np.sum(dist < 2.8) == 0: 
                score = np.sum((dist > 2.8) & (dist < 6.5)) 
                
                if score > max_voisins:
                    max_voisins = score
                    meilleure_poche = point
                    
        return meilleure_poche.tolist()

    except Exception as e:
        print(f"   [PocketFinderV13] Erreur technique : {e}")
        return None

if __name__ == "__main__":
    proteine_test = "1V4S" 
    print(f"--- TEST POCKETFINDER V13 ---")
    print(f"Dossier cible : {DIR_PROTEINES}")
    resultat = trouver_poche_profonde(proteine_test)
    if resultat:
        print(f"Resultat : {resultat}")
    else:
        print(f"Echec du test pour {proteine_test}")
