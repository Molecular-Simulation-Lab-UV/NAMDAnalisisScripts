import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import distances
import multiprocessing as mp
import time

# =========================
# CONFIG
# =========================
topology = "BvPIP2.2.hmr.psf"
trajectory = "/home/dbasaez/CHUPETE/AQPs_RESPALDO/BvPIP2.2/MD1/MD1_BvPIP2.2.tension.gather.dcd"
cutoff = 5.0
n_cores = 6

# =========================
# FUNCIÓN PARALELA
# =========================
def process_frames(args):
    pair, frame_indices = args

    u = mda.Universe(topology, trajectory)

    monomers = {
        "A": u.select_atoms("protein and segid PROA"),
        "B": u.select_atoms("protein and segid PROB"),
        "C": u.select_atoms("protein and segid PROC"),
        "D": u.select_atoms("protein and segid PROD"),
    }

    mono1 = monomers[pair[0]]
    mono2 = monomers[pair[1]]

    res1 = mono1.residues
    res2 = mono2.residues

    # 🔥 MAPEO GLOBAL → LOCAL
    resindex_to_local1 = {res.resindex: i for i, res in enumerate(res1)}
    resindex_to_local2 = {res.resindex: i for i, res in enumerate(res2)}

    atom_to_res1 = np.array([resindex_to_local1[idx] for idx in mono1.resindices])
    atom_to_res2 = np.array([resindex_to_local2[idx] for idx in mono2.resindices])

    contact_map = np.zeros((len(res1), len(res2)))

    for frame in frame_indices:
        u.trajectory[frame]

        # 🔥 DISTANCIA VECTORIAL
        d = distances.distance_array(mono1.positions, mono2.positions)

        contacts = d < cutoff
        idx1, idx2 = np.where(contacts)

        res_i = atom_to_res1[idx1]
        res_j = atom_to_res2[idx2]

        # evitar duplicados por frame
        unique_pairs = set(zip(res_i, res_j))

        for i, j in unique_pairs:
            contact_map[i, j] += 1

    return contact_map


# =========================
# FUNCIÓN PRINCIPAL
# =========================
def compute_contacts_parallel(pair):

    print(f"\nAnalizando {pair[0]}-{pair[1]}...")

    u = mda.Universe(topology, trajectory)
    n_frames = len(u.trajectory)

    frames = np.arange(n_frames)
    chunks = np.array_split(frames, n_cores)

    args = [(pair, chunk) for chunk in chunks]

    with mp.Pool(n_cores) as pool:
        results = pool.map(process_frames, args)

    final_map = np.sum(results, axis=0)
    final_map /= n_frames

    mono1 = u.select_atoms(f"protein and segid PRO{pair[0]}")
    mono2 = u.select_atoms(f"protein and segid PRO{pair[1]}")

    return final_map, mono1.residues, mono2.residues


# =========================
# PARES
# =========================
pairs = [("A","B"), ("A","C"), ("B","D"), ("C","D")]

# =========================
# ⏱️ TIEMPO TOTAL
# =========================
start_time = time.time()

# =========================
# LOOP PRINCIPAL
# =========================
for p in pairs:

    pair_start = time.time()

    cmap, r1, r2 = compute_contacts_parallel(p)

    # ======================================
    # HEATMAP
    # ======================================
    plt.figure(figsize=(6,5))
    plt.imshow(cmap, origin='lower', aspect='auto')
    plt.colorbar(label="Frecuencia de contacto")
    plt.title(f"Contactos {p[0]}-{p[1]}")
    plt.xlabel(f"Residuos {p[1]}")
    plt.ylabel(f"Residuos {p[0]}")
    plt.tight_layout()
    plt.savefig(f"heatmap_{p[0]}_{p[1]}.png", dpi=300)
    plt.close()

    # ======================================
    # MATRIZ CRUDA
    # ======================================
    np.save(f"cmap_{p[0]}_{p[1]}.npy", cmap)

    # ======================================
    # TABLA COMPLETA
    # ======================================
    with open(f"contacts_full_{p[0]}_{p[1]}.dat", "w") as f:
        f.write("i\tj\tRes1\tRes2\tFreq\n")

        for i in range(cmap.shape[0]):
            for j in range(cmap.shape[1]):
                res1 = f"{r1[i].resname}{r1[i].resid}"
                res2 = f"{r2[j].resname}{r2[j].resid}"
                f.write(f"{i}\t{j}\t{res1}\t{res2}\t{cmap[i,j]:.6f}\n")

    # ======================================
    # TOP CONTACTOS
    # ======================================
    threshold = 0.3
    contacts = []

    for i in range(cmap.shape[0]):
        for j in range(cmap.shape[1]):
            if cmap[i,j] > threshold:
                contacts.append((r1[i], r2[j], cmap[i,j]))

    contacts = sorted(contacts, key=lambda x: x[2], reverse=True)

    print(f"\nTop contactos {p[0]}-{p[1]}:")
    for c in contacts[:20]:
        print(f"{c[0].resname}{c[0].resid} - {c[1].resname}{c[1].resid} : {c[2]:.2f}")

    pair_end = time.time()
    print(f"⏱️ Tiempo {p[0]}-{p[1]}: {pair_end - pair_start:.2f} s")


# =========================
# ⏱️ TIEMPO FINAL
# =========================
end_time = time.time()
elapsed = end_time - start_time

print("\n==============================")
print(f"Tiempo total: {elapsed:.2f} segundos")

horas = int(elapsed // 3600)
minutos = int((elapsed % 3600) // 60)
segundos = elapsed % 60

print(f"Tiempo total formateado: {horas}h {minutos}m {segundos:.2f}s")
print("==============================")
