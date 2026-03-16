import vcf
import os
from datetime import datetime

# --- BASE DE CONNAISSANCES PHARMGKB ---
RULES = {
    'rs1799884': {
        'gene': 'GCK',
        'label': 'Diabète Monogénique (MODY 2)',
        'impact': 'Défaut du capteur de glucose pancréatique.',
        'drug_responses': {
            'Metformine': {'score': 20, 'note': 'Peu efficace : glycémie stable mais élevée.'}
        }
    },
    'rs11212617': {
        'gene': 'ATM',
        'label': 'Répondeur Élevé (Type 2)',
        'impact': 'Variante associée à une meilleure réponse à la metformine.',
        'drug_responses': {
            'Metformine': {'score': 90, 'note': 'Excellente réponse attendue.'}
        }
    },
    'rs4124874': {
        'gene': 'SLC22A1',
        'label': 'Transporteur réduit',
        'impact': 'Diminution de l\'absorption de la metformine par le foie.',
        'drug_responses': {
            'Metformine': {'score': 40, 'note': 'Efficacité réduite.'}
        }
    }
}

def analyser_vcf(vcf_path):
    if not os.path.exists(vcf_path):
        print(f"Fichier {vcf_path} introuvable.")
        return
    
    reader = vcf.Reader(filename=vcf_path)
    found_variants = [RULES[record.ID] for record in reader if record.ID in RULES]
    
    # Générer le texte du rapport
    rapport_texte = generer_rapport_string(found_variants)
    
    # Afficher à l'écran
    print(rapport_texte)
    
    # Sauvegarder dans un fichier .txt
    nom_fichier = f"rapport_diabete_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    with open(nom_fichier, "w", encoding="utf-8") as f:
        f.write(rapport_texte)
    print(f"Fichier sauvegardé sous : {nom_fichier}")

def generer_rapport_string(variants):
    lines = []
    lines.append("="*60)
    lines.append("       RAPPORT DE PHARMACOGÉNOMIQUE : DIABÈTE")
    lines.append("="*60)
    lines.append(f"Date de l'analyse : {datetime.now().strftime('%d/%m/%Y %H:%M')}")

    if not variants:
        lines.append("\nAucun variant pertinent détecté.")
    else:
        for var in variants:
            lines.append(f"\nGÈNE : {var['gene']} | DIAGNOSTIC : {var['label']}")
            lines.append(f"IMPACT : {var['impact']}")
            lines.append("-" * 30)
            for drug, data in var['drug_responses'].items():
                lines.append(f" > {drug} : {data['score']}/100 - {data['note']}")
    
    lines.append("\n" + "="*60)
    lines.append("Note : Ce rapport est généré par ordinateur (Projet Diabete).")
    lines.append("="*60)
    return "\n".join(lines)

if __name__ == "__main__":
    analyser_vcf("patient_complet.vcf")
