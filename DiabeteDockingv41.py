import vcf, os, subprocess, warnings, datetime
from fpdf import FPDF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from pocket_finder import trouver_poche_profonde 

warnings.filterwarnings("ignore")
VINA_EXE = "vina.exe" 

# --- BASE DE CONNAISSANCES ---
TABLEAU_STRATEGIQUE = {
    'rs1801278': {'gene': 'INSR', 'pdb': '1IRK', 'panne': 'Binding',    'calcul_vina': True,  'drug': 'Insulin_Mimetic', 'smiles': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C2=O)C=CC=C3O)O)O'},
    'rs5219':    {'gene': 'KCNJ11','pdb': '6C3O', 'panne': 'Binding',    'calcul_vina': True,  'drug': 'Gliclazide',      'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N2CCCC3C2CCCC3'},
    'rs1799884': {'gene': 'GCK',   'pdb': '1V4S', 'panne': 'Activité',   'calcul_vina': False, 'drug': 'Activateur_GCK',  'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N'},
    'rs11212617':{'gene': 'ATM',   'pdb': '7SAY', 'panne': 'Activité',   'calcul_vina': False, 'drug': 'Metformine',      'smiles': 'CN(C)C(=N)N=C(N)N'},
    'rs7903146': {'gene': 'TCF7L2','pdb': '2G4B', 'panne': 'Expression', 'calcul_vina': False, 'drug': 'GLP1_Analogue',   'smiles': 'C1=CC=C(C=C1)CC(C(=O)O)N'},
    'rs4124874': {'gene': 'SLC22A1','pdb': '8H66', 'panne': 'Expression', 'calcul_vina': False, 'drug': 'Bypass_Molecule', 'smiles': 'CN1C2=C(C(=O)N(C1=O)C)N=CN2'}
}

# --- CLASSE PDF ---
class PDFReport(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 14)
        self.cell(0, 10, 'DIABETE-GENE : RAPPORT MEDICAL DE PRECISION', 0, 1, 'C')
        self.ln(5)

# --- FONCTIONS TECHNIQUES ---
def lancer_docking_final(rsid, info):
    ligand = f"{info['drug']}.pdbqt"
    # Préparation ligand (simplifiée pour le flux)
    mol = Chem.MolFromSmiles(info['smiles'])
    mol = Chem.AddHs(mol); AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    prepper = MoleculePreparation(); prepper.prepare(mol)
    with open(ligand, "w") as f: f.write(prepper.write_pdbqt_string())
    
    coords = trouver_poche_profonde(info['pdb'])
    if not coords: return "Erreur Géo"
    cx, cy, cz = coords
    
    conf = f"conf_{rsid}.txt"
    with open(conf, "w") as f:
        f.write(f"receptor = {info['pdb']}.pdbqt\nligand = {ligand}\n")
        f.write(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n")
        f.write("size_x = 18.0\nsize_y = 18.0\nsize_z = 18.0\nexhaustiveness = 8\n")

    try:
        res = subprocess.run([VINA_EXE, "--config", conf], capture_output=True, text=True)
        for line in res.stdout.split('\n'):
            if "   1 " in line: return line.split()[1]
    except: return "N/A"
    return "N/A"

def main():
    print("="*80)
    print("  DIABETE-GENE V4.1 : ANALYSE & GENERATION DU RAPPORT PDF")
    print("="*80)
    
    reader = vcf.Reader(filename="patient_expert.vcf")
    resultats_pour_pdf = {}

    for record in reader:
        if record.ID in TABLEAU_STRATEGIQUE:
            info = TABLEAU_STRATEGIQUE[record.ID]
            print(f"\n[DIAGNOSTIC] {record.ID} ({info['gene']})")
            
            score = "N/A"
            if info['calcul_vina']:
                score = lancer_docking_final(record.ID, info)
                print(f"   >>> Score Vina : {score} kcal/mol")
            else:
                print(f"   >>> Panne {info['panne']} (Calcul ignoré)")
            
            # Stockage pour le PDF
            info['score'] = score
            resultats_pour_pdf[record.ID] = info

    # --- GENERATION PDF ---
    pdf = PDFReport()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Synthese des analyses :", 0, 1)
    pdf.set_font("Arial", '', 10)
    
    for rsid, res in resultats_pour_pdf.items():
        txt = f"Mutation: {rsid} | Gene: {res['gene']} | Panne: {res['panne']} | Medicament: {res['drug']}"
        if res['calcul_vina']: txt += f" | Affinite: {res['score']} kcal/mol"
        pdf.multi_cell(0, 8, txt, 1)
        pdf.ln(2)

    pdf.output("Rapport_Final_Patient.pdf")
    print("\n" + "="*80)
    print("✅ ANALYSE TERMINEE : Ouvrez 'Rapport_Final_Patient.pdf'")

if __name__ == "__main__":
    main()
