import vcf
import os
from datetime import datetime

# --- BASE PHARMGKB ENRICHIE ---
# Scores : 0-100 (Efficacité/Sécurité)
RULES = {
    # --- DIAGNOSTIC & RÉPONSE MÉDICAMENTS ---
    'rs1799884': {
        'gene': 'GCK',
        'label': 'MODY 2 (Glucokinase)',
        'impact': 'Seuil de détection du glucose élevé.',
        'drugs': {
            'Metformine': {'score': 15, 'note': 'Peu d\'effet sur la glycémie stable de MODY2.'},
            'Sulfonylurées': {'score': 10, 'note': 'Risque inutile d\'hypoglycémie.'},
            'Hygiène de vie': {'score': 95, 'note': 'Traitement de première intention recommandé.'}
        }
    },
    'rs5219': {
        'gene': 'KCNJ11',
        'label': 'Diabète Néonatal / MODY',
        'impact': 'Mutation des canaux potassiques du pancréas.',
        'drugs': {
            'Sulfonylurées': {'score': 95, 'note': 'Réponse exceptionnelle ! Souvent supérieur à l\'insuline.'},
            'Metformine': {'score': 40, 'note': 'Efficacité modérée.'}
        }
    },
    'rs7903146': {
        'gene': 'TCF7L2',
        'label': 'Prédisposition Type 2',
        'impact': 'Altération de la sécrétion d\'insuline et effet incrétine.',
        'drugs': {
            'GLP-1 (Ozempic/Victoza)': {'score': 85, 'note': 'Très bonne réponse aux agonistes GLP-1.'},
            'Sulfonylurées': {'score': 50, 'note': 'Efficacité réduite chez les porteurs du variant T.'}
        }
    },
    # --- PHARMACOGÉNOMIQUE PURE (TRANSPORT & MÉTABOLISME) ---
    'rs10305420': {
        'gene': 'SLC22A1',
        'label': 'Transporteur OCT1 réduit',
        'impact': 'Diminution du transport de la metformine vers le foie.',
        'drugs': {
            'Metformine': {'score': 35, 'note': 'Réponse suboptimale (Niveau PharmGKB 1A).'}
        }
    },
    'rs12248560': {
        'gene': 'CYP2C19',
        'label': 'Métaboliseur Rapide',
        'impact': 'Dégradation accélérée de certains antidiabétiques oraux.',
        'drugs': {
            'Gliclazide': {'score': 60, 'note': 'Élimination rapide, peut nécessiter un ajustement de dose.'}
        }
    }
}

def analyser_patient_expert(vcf_path):
    if not os.path.exists(vcf_path):
        print(f"Fichier {vcf_path} absent.")
        return

    reader = vcf.Reader(filename=vcf_path)
    found_variants = [RULES[record.ID] for record in reader if record.ID in RULES]
    
    # Construction du rapport
    output = []
    output.append("="*70)
    output.append("   SYSTÈME EXPERT DE DIAGNOSTIC GÉNÉTIQUE & PHARMACOLOGIQUE")
    output.append("="*70)
    output.append(f"Analyse générée le : {datetime.now().strftime('%d/%m/%Y %H:%M')}")

    if not found_variants:
        output.append("\nAucun marqueur génétique majeur détecté dans les zones ciblées.")
    else:
        for var in found_variants:
            output.append(f"\n[GÈNE : {var['gene']}] -> {var['label']}")
            output.append(f"DESCRIPTION : {var['impact']}")
            output.append("-" * 40)
            output.append("SCORES D'EFFICACITÉ PERSONNALISÉS :")
            
            for drug, data in var['drugs'].items():
                status = "PRESCRIRE" if data['score'] >= 80 else "PRÉCAUTION" if data['score'] >= 40 else "ÉVITER"
                output.append(f"  * {drug:15} : {data['score']}/100 [{status}]")
                output.append(f"    Note : {data['note']}")

    output.append("\n" + "="*70)
    output.append("AVERTISSEMENT : Ce score est basé sur les données PharmGKB.")
    output.append("Consultez toujours un endocrinologue avant toute modification.")
    output.append("="*70)

    # Exportation
    resultat_final = "\n".join(output)
    print(resultat_final)
    
    filename = f"expert_report_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
    with open(filename, "w", encoding="utf-8") as f:
        f.write(resultat_final)
    print(f"\n[SUCCESS] Rapport exporté : {filename}")

if __name__ == "__main__":
    analyser_patient_expert("patient_complet.vcf")
