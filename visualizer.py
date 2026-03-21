import matplotlib.pyplot as plt
import pandas as pd
import os
from joblib import load
from sklearn.metrics import ConfusionMatrixDisplay, roc_curve, auc

def produce_figures():
    if not os.path.exists("Resultats/trained_model.joblib"):
        print("❌ Erreur : Modèle introuvable.")
        return

    model = load("Resultats/trained_model.joblib")
    test_data = pd.read_csv("Resultats/data_test_real.csv")
    
    X_test = test_data.drop('label', axis=1)
    y_test = test_data['label']

    print("📊 Génération des graphiques...")

    # --- 1. MATRICE DE CONFUSION ---
    plt.figure(figsize=(6, 5))
    # On force les labels [0, 1] pour éviter l'erreur de shape
    ConfusionMatrixDisplay.from_predictions(y_test, model.predict(X_test), 
                                            display_labels=['Sensible', 'Résistant'],
                                            labels=[0, 1], cmap='Blues')
    plt.title("Matrice de Confusion (Validation)")
    plt.savefig("Resultats/Figure_1_CM.tiff", dpi=300)
    plt.close()

    # --- 2. COURBE ROC (SÉCURISÉE) ---
    plt.figure(figsize=(7, 6))
    
    # Vérification si on a bien les deux classes (0 et 1)
    if len(model.classes_) > 1:
        probs = model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, probs)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, color='darkorange', label=f"AUC = {roc_auc:.2f}")
        plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
        plt.title("Courbe ROC - Performance du Modèle")
    else:
        # Cas où le modèle n'a appris qu'une seule classe
        plt.text(0.5, 0.5, "Données insuffisantes pour ROC\n(Une seule classe détectée)", 
                 ha='center', va='center', fontsize=12, color='red')
        plt.title("Courbe ROC (Non disponible)")
        print("⚠️ Note : Une seule classe détectée dans les données. La courbe ROC est désactivée.")

    plt.legend(loc="lower right")
    plt.savefig("Resultats/Figure_2_ROC.tiff", dpi=300)
    plt.close()
    
    print("✅ Figures TIFF (300 DPI) générées dans /Resultats.")
