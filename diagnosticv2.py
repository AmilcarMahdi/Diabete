import vcf
import os

# --- BASE DE CONNAISSANCES PHARMGKB (Simplifiée pour le Diabète) ---
RULES = {
    # Diagnostic de Type / Risque
    'rs1799884': {
        'gene': 'GCK',
        'label': 'Diabète Monogénique (MODY 2)',
        'impact': 'Défaut du capteur de glucose pancréatique.',
        'drug_responses': {
            'Metformine': {'score': 20, 'note': 'Peu efficace : glycémie stable mais élevée.'},
            'Sulfonylurées': {'score': 10, 'note': 'Non recommandé (risque d\'hypoglycémie).'}
        }
    },
    # Pharmacogénomique : Réponse à la Metformine (Gène ATM)
    'rs11212617': {
        'gene': 'ATM',
        'label': 'Répondeur Élevé (Type 2)',
        'impact': 'Variante associée à une meilleure réponse glycémique à la metformine.',
        'drug_responses': {
            'Metformine': {'score': 90, 'note': 'Excellente réponse attendue.'}
        }
    },
    # Pharmacogénomique : Transporteur SLC22A1 (OCT1)
    'rs4124874': {
        'gene': 'SLC22A1',
        'label': 'Transporteur réduit',
        'impact': 'Diminution de l\'absorption de la metformine par le foie.',
        'drug_responses': {
            'Metformine': {'score': 40, 'note': 'Efficacité réduite : nécessite souvent une dose plus forte ou alternative.'}
        }
    }
}

def analyser_vcf(vcf_path):
    if not os.path.exists(vcf_path):
        print(f"Fichier {vcf_path} introuvable.")
        return

    reader = vcf.Reader(filename=vcf_path)
    found_variants = []

    for record in reader:
        if record.ID in RULES:
            found_variants.append(RULES[record.ID])

    generer_rapport(found_variants)

def generer_rapport(variants):
    print("\n" + "="*60)
    print("       RAPPORT DE PHARMACOGÉNOMIQUE : DIABÈTE")
    print("="*60)

    if not variants:
        print("Aucun variant pathogène ou pharmacogénomique détecté.")
        return

    for var in variants:
        print(f"\nGENE : {var['gene']} | DIAGNOSTIC : {var['label']}")
        print(f"IMPACT BIOLOGIQUE : {var['impact']}")
        print("-" * 30)
        print("PRÉVISIONS DE RÉPONSES AUX MÉDICAMENTS :")
        
        for drug, data in var['drug_responses'].items():
            color = "OK" if data['score'] >= 70 else "ATTENTION"
            print(f" > [{color}] {drug} : {data['score']}/100")
            print(f"   Note : {data['note']}")
    
    print("\n" + "="*60)
    print("Note : Ce rapport est généré par IA à des fins de recherche.")
    print("="*60 + "\n")

if __name__ == "__main__":
    analyser_vcf("patient_complet.vcf")
