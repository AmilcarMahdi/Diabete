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
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    
    atoms = np.array(coords)
    # Définition des limites de la protéine
    min_c, max_c = np.min(atoms, axis=0), np.max(atoms, axis=0)
    
    # Génération de points tests (Grid Search)
    grid = np.random.uniform(min_c, max_c, (2000, 3))
    
    meilleure_poche = np.mean(atoms, axis=0)
    score_max = 0
    
    for point in grid:
        dist = np.linalg.norm(atoms - point, axis=1)
        # Un point dans une poche est vide à 3A mais entouré à 10A
        if np.sum(dist < 3.0) == 0: 
            densite_entourage = np.sum((dist > 3.5) & (dist < 12.0))
            if densite_entourage > score_max:
                score_max = densite_entourage
                meilleure_poche = point
                
    return meilleure_poche.tolist()
    
# --- BLOC DE TEST À AJOUTER ---
if __name__ == "__main__":
    # Test sur une de tes protéines (ex: 1V4S)
    proteine_test = "1V4S" 
    print(f"Recherche de la serrure pour {proteine_test}...")
    
    resultat = trouver_poche_profonde(proteine_test)
    
    if resultat:
        print(f"✅ Serrure trouvée ! Coordonnées : {resultat}")
    else:
        print("❌ Erreur : Fichier .pdb introuvable dans le dossier.")
