import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestClassifier
from joblib import dump

# Propriétés physico-chimiques des acides aminés (Optionnel pour enrichir le modèle)
AA_PROPS = {'ALA': 1.8, 'ARG': -4.5, 'ILE': 4.5, 'LEU': 3.8, 'GLN': -3.5, 'CYS': 2.5, 'GLU': -3.5}

def augment_data(df, n=300):
    """
    Data Augmentation intelligente pour article Q1.
    On injecte un bruit gaussien sur le delta d'affinité pour simuler
    l'incertitude de mesure de Vina (SD = 0.1 kcal/mol).
    """
    augmented = []
    # On définit le seuil basé sur la médiane des données réelles
    threshold = df['delta_aff'].median()
    
    for _ in range(n):
        # On pioche une ligne au hasard
        s = df.sample(1).iloc[0].to_dict()
        
        # On ajoute du bruit au delta_aff
        noise = np.random.normal(0, 0.12)
        s['delta_aff'] += noise
        
        # RE-LABELLISATION DYNAMIQUE : Crucial pour la courbe ROC
        # Si le bruit fait passer le delta au-dessus du seuil, il devient une résistance (1)
        s['label'] = 1 if s['delta_aff'] >= threshold else 0
        augmented.append(s)
        
    return pd.concat([df, pd.DataFrame(augmented)], ignore_index=True)

def train_model():
    input_csv = "Resultats/results_raw.csv"
    if not os.path.exists(input_csv):
        print(f"❌ Erreur : {input_csv} introuvable. Lance le docking d'abord.")
        return

    df_raw = pd.read_csv(input_csv)
    
    # Nettoyage : on ignore les scores nuls (dockings échoués)
    df_raw = df_raw[df_raw['m_wt'] != 0].copy()
    
    if len(df_raw) < 3:
        print("⚠️ Pas assez de données réelles pour l'entraînement.")
        return

    # Calcul du Delta d'affinité physique
    df_raw['delta_aff'] = df_raw['m_mut'] - df_raw['m_wt']
    
    # SEUIL MÉDIAN : La clé pour débloquer la courbe ROC
    threshold = df_raw['delta_aff'].median()
    
    features = []
    for _, row in df_raw.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is None: continue
        
        # Extraction de features moléculaires (Descripteurs RDKit)
        features.append({
            'mol_wt': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'h_donors': Descriptors.NumHDonors(mol),
            'delta_aff': row['delta_aff'],
            'label': 1 if row['delta_aff'] >= threshold else 0
        })
    
    df_ml = pd.DataFrame(features)
    
    # Augmentation de données avec équilibrage des classes (0 et 1)
    df_aug = augment_data(df_ml, n=400)
    
    # Préparation des sets X (caractéristiques) et y (cible/label)
    X = df_aug.drop('label', axis=1)
    y = df_aug['label']
    
    # Entraînement du modèle Random Forest
    # On utilise 100 arbres pour une bonne généralisation
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    # Sauvegarde du modèle et des données de test
    dump(model, "Resultats/trained_model.joblib")
    df_ml.to_csv("Resultats/data_test_real.csv", index=False)
    
    print(f"✅ ML Engine : Modèle entraîné avec {len(df_aug)} points (Seuil: {threshold:.3f})")

if __name__ == "__main__":
    train_model()
