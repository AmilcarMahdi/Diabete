import numpy as np
import os

def trouver_poche_profonde(pdb_id):
    # FORCE LE RÉSULTAT À ÊTRE TOUJOURS LE MÊME
    np.random.seed(42) 
    
    pdb_file = f"{pdb_id}.pdb"
    if not os.path.exists(pdb_file):
        return None
        
    coords = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    
    atoms = np.array(coords)
    min_c, max_c = np.min(atoms, axis=0), np.max(atoms, axis=0)
    
    # GÉNÈRE PLUS DE POINTS POUR LA PRÉCISION
    grid = np.random.uniform(min_c, max_c, (5000, 3)) 
    
    # --- LA CORRECTION EST ICI : INITIALISATION ---
    meilleure_poche = np.mean(atoms, axis=0)
    max_voisins = 0 # <--- INDISPENSABLE
    
    for point in grid:
        dist = np.linalg.norm(atoms - point, axis=1)
        # Condition de vide : pas d'atome à moins de 2.5 Angströms
        if np.sum(dist < 2.5) == 0: 
            # Score de densité : atomes entre 2.5 et 6.0 Angströms
            score = np.sum((dist > 2.5) & (dist < 6.0)) 
            if score > max_voisins:
                max_voisins = score
                meilleure_poche = point
                
    return meilleure_poche.tolist()

# --- BLOC DE TEST ---
if __name__ == "__main__":
    proteine_test = "1V4S" 
    print(f"--- RECHERCHE DE LA SERRURE PROFONDE ({proteine_test}) ---")
    coords = trouver_poche_profonde(proteine_test)
    if coords:
        print(f"✅ FOND DE SERRURE TROUVÉ ! Coordonnées : {coords}")
    else:
        print("❌ Erreur : Fichier PDB manquant.")
