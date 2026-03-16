import vcf  # On utilise PyVCF3 installé précédemment
import pandas as pd

# --- ÉTAPE 1 : Base de connaissances simplifiée (PharmGKB style) ---
# On définit l'impact des gènes sur le diagnostic et le traitement
CONNAISSANCES = {
    'GCK': {'type': 'MODY2', 'metformine_eff': 0.5, 'description': 'Défaut de détection du glucose'},
    'HNF1A': {'type': 'MODY3', 'metformine_eff': 0.8, 'description': 'Sensibilité aux sulfonylurées'},
    'TCF7L2': {'type': 'Type 2', 'metformine_eff': 1.2, 'description': 'Risque élevé de diabète de type 2'},
    'SLC22A1': {'type': 'Pharmacogène', 'metformine_eff': 0.4, 'description': 'Transporteur de metformine réduit'}
}

def analyser_patient(vcf_file):
    print(f"--- Analyse du fichier : {vcf_file} ---")
    
    # Simulation : Dans un vrai VCF, on chercherait les variants spécifiques.
    # Ici, on simule la détection de deux gènes mutés chez le patient.
    genes_detectes = ['TCF7L2', 'SLC22A1'] 
    
    score_metformine = 100 # Score de base
    diagnostics = []

    for gene in genes_detectes:
        if gene in CONNAISSANCES:
            info = CONNAISSANCES[gene]
            diagnostics.append(f"Mutation trouvée sur {gene} ({info['type']}) : {info['description']}")
            
            # Ajustement du score de médicament (multiplication)
            score_metformine *= info['metformine_eff']

    return diagnostics, score_metformine

# --- ÉTAPE 2 : Exécution ---
resultats, score_final = analyser_patient("mon_genome.vcf")

print("\n[Rapport de Diagnostic]")
for res in resultats:
    print(f"- {res}")

print(f"\n[Score Médicament Personnalisé]")
print(f"Efficacité estimée de la Metformine : {score_final:.1f}/100")

if score_final < 50:
    print("Conseil : Ce patient pourrait nécessiter une alternative à la Metformine.")
