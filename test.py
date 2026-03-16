import subprocess

# On teste sur la première mutation qui a échoué
cmd = ["vina.exe", "--config", "conf_rs1799884.txt", "--out", "test_out.pdbqt", "--log", "test_log.txt"]

print("--- TEST DEBOGAGE VINA ---")
result = subprocess.run(cmd, capture_output=True, text=True)

print("STDOUT (Sortie classique) :", result.stdout)
print("\nSTDERR (L'ERREUR REELLE) :", result.stderr)
