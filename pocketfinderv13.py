import numpy as np
import os

DIR_PROTEINES = "Proteines"

def trouver_poche_profonde(pdb_input):
    # 1. Gestion du chemin
    if os.path.exists(pdb_input):
        pdb_path = pdb_input
    else:
        pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_input}.pdbqt")
        if not os.path.exists(pdb_path):
            pdb_path = os.path.join(DIR_PROTEINES, f"{pdb_input}.pdb")

    if not os.path.exists(pdb_path): return None
        
    coords = []
    # 2. Extraction des coordonnées ATOM uniquement (on ignore le vide)
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: continue
    
    if not coords: return None
    atoms = np.array(coords)

    # 3. ALGORITHME DE DENSITÉ (Recherche du point le plus "entouré")
    # On cherche l'atome qui a le plus de voisins dans un rayon de 10A
    best_atom_idx = 0
    max_neighbors = 0
    
    # Pour accélérer, on échantillonne un atome sur 5
    for i in range(0, len(atoms), 5):
        # Distance entre l'atome i et tous les autres
        diff = atoms - atoms[i]
        dist_sq = np.sum(diff**2, axis=1)
        # On compte combien d'atomes sont entre 3.5A et 10A (la "coquille" de la poche)
        neighbors = np.sum((dist_sq > 12.25) & (dist_sq < 100.0))
        
        if neighbors > max_neighbors:
            max_neighbors = neighbors
            best_atom_idx = i
            
    # Le centre de la poche est l'atome le plus enfoui
    poche_coords = atoms[best_atom_idx].tolist()
    
    return poche_coords

if __name__ == "__main__":
    test_id = "1V4S"
    print(f"--- TEST POCKETFINDER V15 (DENSITÉ) ---")
    res = trouver_poche_profonde(test_id)
    print(f"Coordonnées trouvées : {res}")
