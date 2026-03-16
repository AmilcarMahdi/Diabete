import numpy as np
import os

def trouver_poche_profonde(pdb_id):
    pdb_file = f"{pdb_id}.pdb"
    if not os.path.exists(pdb_file):
        return None
        
    coords = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # Extraction des coordonnées X, Y, Z du format PDB
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    
    atoms = np.array(coords)
    # Définition des limites de la protéine pour la zone de recherche
    min_c, max_c = np.min(atoms, axis=0), np.max(atoms, axis=0)
    
    # Génération de 3000 points tests (on augmente la densité pour trouver le fond)
    grid = np.random.uniform(min_c, max_c, (3000, 3))
    
    meilleure_poche = np.mean(atoms, axis=0)
    score_max = 0
    
    for point in grid:
        dist = np.linalg.norm(atoms - point, axis=1)
        
        # CONDITION 1 : Le point doit être dans un VIDE (pas de collision avec un atome)
        # On passe à 2.8A pour s'assurer que même un petit ligand rentre
        if np.sum(dist < 2.8) == 0: 
            
            # CONDITION 2 : Densité d'entourage serrée (Le "Fond du trou")
            # En réduisant le rayon à 7.0A, on pénalise les entrées larges 
            # et on favorise les cavités étroites et profondes.
            voisins_fond = np.sum((dist > 2.8) & (dist < 7.0))
            
            if voisins_fond > score_max:
                score_max = voisins_fond
                meilleure_poche = point
                
    return meilleure_poche.tolist()

# --- BLOC DE TEST POUR AFFICHAGE ---
if __name__ == "__main__":
    proteine_test = "1V4S" 
    print(f"--- RECHERCHE DE LA SERRURE PROFONDE ({proteine_test}) ---")
    
    coords_finales = trouver_poche_profonde(proteine_test)
    
    if coords_finales:
        print(f"✅ FOND DE SERRURE TROUVÉ !")
        print(f"Coordonnées : {coords_finales}")
        print("\nCopie ces chiffres dans PyMOL pour vérifier :")
        print(f"pseudoatom ma_serrure, pos={coords_finales}")
    else:
        print(f"❌ Erreur : {proteine_test}.pdb est absent du dossier.")
