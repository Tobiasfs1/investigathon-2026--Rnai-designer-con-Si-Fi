# -*- coding: utf-8 -*-
# Python 2. Guardar en la misma carpeta que general_helpers.py
#
# Uso:
#   python region_verde.py salida_sifi.json mi_gen.fasta 21
#
# (json de si-Fi, fasta del gen, tamano de siRNA que elegiste)

import sys
import general_helpers

if len(sys.argv) != 4:
    print("Uso: python region_verde.py salida_sifi.json mi_gen.fasta 21")
    sys.exit(1)

f_in = sys.argv[1]
fasta = sys.argv[2]
sirna_size = int(sys.argv[3])

seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()

off, main, eff, main_hits = general_helpers.get_target_data(f_in, sirna_size)

verde = set()
for v in main.values():
    verde |= v

if not verde:
    print("No hay region verde (main target sin hits).")
    sys.exit(1)

ini = min(verde)   # 0-based, como el eje X del grafico
fin = max(verde)
largo = fin - ini + 1

# posiciones cubiertas por siRNAs eficientes dentro de la region
n_eff = sum(1 for v in eff.values() for p in v if ini <= p <= fin)

print("Region verde: %d-%d (1-based, inclusivo)" % (ini + 1, fin + 1))
print("Largo: %d nt" % largo)
print("Cobertura eficiente: %d posiciones (%.2f por nt)" % (n_eff, n_eff / float(largo)))
print("")
print(seq[ini:fin + 1])

# FASTA listo para el filtro de off-targets
out = "region_verde_%d-%d.fasta" % (ini + 1, fin + 1)
with open(out, "w") as fh:
    fh.write(">region_verde_%d-%d\n%s\n" % (ini + 1, fin + 1, seq[ini:fin + 1]))
print("")
print("Guardado en: " + out)
