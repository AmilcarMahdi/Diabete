import os, subprocess, re, numpy as np, pandas as pd
from Bio.PDB import PDBParser, PDBIO
from rdkit import Chem               # <-- AJOUTE CETTE LIGNE
from rdkit.Chem import AllChem      # <-- AJOUTE CETTE LIGNE
from pocketfinderv13 import trouver_poche_profonde


# --- CONFIGURATION ---
VINA_EXE = "vina.exe"
DIRS = ["Proteines", "Ligands", "Resultats"]
for d in DIRS: os.makedirs(d, exist_ok=True)

# Mapping pour la mutation (Acides Aminés 1 lettre vers 3 lettres)
AA1to3 = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS',
    'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
    'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'
}

def run_multi_docking(receptor_pdb, ligand_pdbqt, center, n_runs=3):
    if not os.path.exists(VINA_EXE): return (0, 0)

    # 1. Chemins ABSOLUS pour Windows
    receptor_abs = os.path.abspath(receptor_pdb)
    ligand_abs = os.path.abspath(ligand_pdbqt)
    receptor_pdbqt = receptor_abs + "qt"
    # 2. GESTION DU CENTRE ET DE LA TAILLE (POCKET FINDER)
    if center is None:
        print(f"   🔍 PocketFinder V13 recherche la poche sur {os.path.basename(receptor_pdb)}...")
        center_coords = trouver_poche_profonde(receptor_pdb)
        if center_coords is None:
            print("   ⚠️ Poche introuvable, utilisation du centre géométrique par défaut.")
            center_coords = [0.0, 0.0, 0.0] 
        size = 40  # Boîte large pour le Blind Docking
    else:
        center_coords = center
        size = 25  # Boîte ciblée si le centre est connu (ex: GCK)

    c_x, c_y, c_z = center_coords
    
    # 2. Conversion PDB -> PDBQT (NETTOYAGE COMPLET)
    with open(receptor_abs, 'r') as f_in, open(receptor_pdbqt, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                element = line[12:14].strip().upper()
                if element in ["ZN", "MG", "FE", "H"]: continue 
                f_out.write(f"{line[:54]:54s}  1.00  0.00                        {element:2s}\n")

    scores = []
    for i in range(n_runs):
        seed = np.random.randint(0, 10000)
        cmd = [
            VINA_EXE, 
            "--receptor", receptor_pdbqt, 
            "--ligand", ligand_abs,
            "--center_x", f"{c_x:.3f}", "--center_y", f"{c_y:.3f}", "--center_z", f"{c_z:.3f}",
            "--size_x", str(size), "--size_y", str(size), "--size_z", str(size), 
            "--exhaustiveness", "8", "--seed", str(seed)
        ]
        
        res = subprocess.run(cmd, capture_output=True, text=True)

        # 3. EXTRACTION ROBUSTE DU SCORE (Correction du 0.00)
        # On nettoie les retours à la ligne Windows \r\n
        output_clean = res.stdout.replace('\r', '')
        
        if "-----+------------+----------+----------" in output_clean:
            # On coupe le texte juste après la ligne de séparation du tableau
            parts = output_clean.split("-----+------------+----------+----------")
            table_rows = parts[-1].strip().split('\n')
            if table_rows:
                # La première ligne après les tirets contient le score du Mode 1
                # On extrait tous les nombres de cette ligne
                row_values = re.findall(r"[-+]?\d+\.\d+", table_rows[0])
                if row_values:
                    # Le premier nombre trouvé est TOUJOURS l'affinité
                    val_score = float(row_values[0])
                    scores.append(val_score)
                    print(f"      -> Run {i+1}: {val_score} kcal/mol")
        else:
            if "mode |   " in output_clean: # Secours si le formatage change légèrement
                found = re.findall(r"(?<!\d)-\d+\.\d+", output_clean)
                if found: scores.append(float(found[0]))

    return (np.mean(scores), np.std(scores)) if scores else (0.0, 0.0)


def generate_raw_results(mutations, db):
    """
    Pipeline complet : Mutation -> Génération Ligand 3D compatible Vina -> Docking REEL -> Analyse Delta.
    """
    results = []
    parser = PDBParser(QUIET=True)
    
    for rsid, mut in mutations.items():
        gene = mut['gene']
        if gene not in db:
            print(f"⏩ Gène {gene} non présent dans la DB. Passage...")
            continue
        
        pdb_id = db[gene]['pdb']
        pdb_wt_path = os.path.join("Proteines", f"{pdb_id}.pdb")
        
        if not os.path.exists(pdb_wt_path):
            print(f"❌ Fichier {pdb_wt_path} introuvable dans /Proteines.")
            continue

        # --- 1. CHARGEMENT ET MUTATION ---
        structure = parser.get_structure(gene, pdb_wt_path)
        model = structure[0]
        
        if mut['chain'] not in model:
            print(f"❌ Chaîne {mut['chain']} introuvable dans {pdb_id}")
            continue
            
        chain = model[mut['chain']]
        
        # Recherche du résidu par numéro (id[1])
        target_residue = None
        for res in chain:
            if res.id[1] == mut['pos']: 
                target_residue = res
                break
        
        if target_residue is None:
            print(f"❌ Résidu {mut['pos']} introuvable dans la chaîne {mut['chain']} de {pdb_id}")
            continue

        # Application de la mutation
        old_resname = target_residue.resname
        new_resname = AA1to3.get(mut['alt'], 'ALA')
        target_residue.resname = new_resname
        
        # Sauvegarde du PDB Mutant
        pdb_mut_path = os.path.join("Proteines", f"{pdb_id}_{rsid}_mut.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_mut_path)
        print(f"✅ Mutant généré : {old_resname}{mut['pos']}{new_resname} ({rsid})")

        # --- 2. BOUCLE DE DOCKING SUR LES VRAIS LIGANDS ---
        center = db[gene].get('center')
        for drug, smiles in db[gene]['drugs'].items():
            lig_path = os.path.join("Ligands", f"{drug}.pdbqt")
            
            # --- GÉNÉRATION ET NETTOYAGE STRICT DU LIGAND (FORMAT ROOT/ENDROOT) ---
            if not os.path.exists(lig_path) or os.path.getsize(lig_path) < 100:
                print(f"   🛠️  Génération 3D du ligand : {drug}...")
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                        pdb_block = Chem.MolToPDBBlock(mol)
                        
                        # Conversion PDB -> PDBQT compatible Vina (Marqueurs ROOT obligatoires)
                        with open(lig_path, "w") as f_lig:
                            f_lig.write("REMARK  LIGAND RIGIDE GENERE\n")
                            f_lig.write("ROOT\n") # Balise indispensable pour Vina
                            for line in pdb_block.splitlines():
                                if line.startswith(("ATOM", "HETATM")):
                                    element = line[12:14].strip().upper()
                                    # On force 'ATOM', on renomme le résidu en 'LIG' et on aligne l'élément en col 77
                                    clean_line = f"ATOM  {line[6:17]}LIG{line[20:54]}  1.00  0.00           {element:2s}\n"
                                    f_lig.write(clean_line)
                            f_lig.write("ENDROOT\n") # Fin du bloc rigide
                            f_lig.write("TORSDOF 0\n") # Degrés de liberté (0 = rigide)
                    else:
                        print(f"   ❌ SMILES invalide pour {drug}")
                        continue
                except Exception as e:
                    print(f"   ❌ Erreur RDKit pour {drug}: {e}")
                    continue

            print(f"   🧪 Docking REEL : {drug} sur {gene}...")

            # Calculs via Vina (WT vs Mutant)
            m_wt, s_wt = run_multi_docking(pdb_wt_path, lig_path, center)
            m_mut, s_mut = run_multi_docking(pdb_mut_path, lig_path, center)
            
            # Stockage des données
            results.append({
                'rsid': rsid, 'gene': gene, 'drug': drug, 'smiles': smiles,
                'm_wt': m_wt, 's_wt': s_wt, 'm_mut': m_mut, 's_mut': s_mut,
                'ref': mut['ref'], 'alt': mut['alt'], 'pos': mut['pos'],
                'chain': mut['chain']
            })
    
    # --- 3. ANALYSE ET LABELLISATION ---
    if not results:
        print("⚠️ Aucune simulation n'a abouti.")
        return

    df = pd.DataFrame(results)
    df['delta_aff'] = df['m_mut'] - df['m_wt']
    
    # Labellisation par la médiane
    #seuil = df['delta_aff'].median()
    df['label'] = (df['delta_aff'] > df['delta_aff'].median()).astype(int)

    
    # --- AFFICHAGE FINAL ---
    print("\n" + "="*85)
    print(f"{'GENE':<8} | {'DRUG':<15} | {'WT (kcal)':<12} | {'MUT (kcal)':<12} | {'DELTA':<8} | {'LABEL'}")
    print("-" * 85)
    for _, r in df.iterrows():
        print(f"{r['gene']:<8} | {r['drug']:<15} | {r['m_wt']:<12.2f} | {r['m_mut']:<12.2f} | {r['delta_aff']:<8.2f} | {int(r['label'])}")
    print("="*85)

    output_path = os.path.join("Resultats", "results_raw.csv")
    df.to_csv(output_path, index=False)
    print(f"\n📁 {len(df)} simulations enregistrées dans {output_path}")


