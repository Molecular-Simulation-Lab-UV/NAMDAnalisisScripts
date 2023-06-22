# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:48:32 2023

@author: AFC
"""
import sys
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
HydrogenBondAnalysis as HBA)
from datetime import datetime
import numpy as np
import pandas as pd

start = datetime.now()

path = sys.argv[1]
print("---------------->", path)
pdb = f'{path}FaPIP21.POPC.Wat.box.ion.pdb'
psf = f'{path}FaPIP21.POPC.Wat.box.ion.psf'
dcd = sys.argv[2]
out = sys.argv[3]

print("---------------->", dcd)
u = MDAnalysis.Universe(psf, dcd)

protres = "51 59 83 87 91 103 104 105 107 110 154 155 157 158 181 185 188 189 190 191 193 197 208 212 216 224 225 226 227 228 231"

hbonds = HBA(universe=u,
             between=['resname TIP3', f'protein and resid {protres}'])
protein_hydrogens_sel = hbonds.guess_hydrogens(f"protein and resid {protres}")
protein_acceptors_sel = hbonds.guess_acceptors(f"protein and resid {protres}") #segid "A" para especificar cadena A

water_hydrogens_sel = "resname TIP3 and name H1 H2"
water_acceptors_sel = "resname TIP3 and name OH2"
hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({water_hydrogens_sel} and around 10 (protein and resid {protres}))" # ELIMINADO --->  and around 10 not resname TIP3
hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({water_acceptors_sel} and around 10 (protein and resid {protres}))"
results = hbonds.run()

elapsed = datetime.now() - start

print("Tiempo de ejecucion de hbonds de MDAnalysis: " ,elapsed, "horas") 

pdbs = dict()
for c in "ABCD":
    pdbs[c] = f'/home/acaviglia/pfs/RefFrames/{c}.ref.pdb'


pdbA = MDAnalysis.Universe(psf, pdbs["A"])
pdbB = MDAnalysis.Universe(psf, pdbs["B"])
pdbC = MDAnalysis.Universe(psf, pdbs["C"])
pdbD = MDAnalysis.Universe(psf, pdbs["D"])

frames = [int(a) for a in results.hbonds[:,0]]
donores = [int(a) for a in results.hbonds[:,1]]
aceptores = [int(a) for a in results.hbonds[:,3]]

datoms = np.array([donores])
datoms = datoms.T
aatoms = np.array([aceptores])
aatoms = aatoms.T
fdata = np.array(frames).T

isA = np.zeros((len(aatoms), 1), dtype=np.int8)
isD = datoms >-1
aar = np.concatenate((aatoms, isA), axis=1)
dar = np.concatenate((datoms, isD), axis=1)
far = np.concatenate((fdata,fdata), axis=0)
adar = np.concatenate((dar, aar), axis=0)

from tqdm import tqdm
zlist = []
clist = []
flist = []
for n,i in tqdm(enumerate(adar[:, 0])):
    i = int(i)
    c = pdbA.atoms[i].segid
    if c == "A":
        atomo = pdbA.select_atoms(f"id {i}")
    elif c == "B":
        atomo = pdbB.select_atoms(f"id {i}")
    elif c == "C":
        atomo = pdbC.select_atoms(f"id {i}")
    elif c == "D":
        atomo = pdbD.select_atoms(f"id {i}")
    else:
        pass
    flist.append(far[n])
    zlist.append(atomo.positions[0,2])
    clist.append(c)

zar = np.array([zlist])
car = np.array([clist])

adar = np.concatenate((np.array([flist]).T, adar, car.T, zar.T), axis=1)

column_names = ["frame", "id", "isD", "chain", "z"]
addf = pd.DataFrame(adar, index=None, columns=column_names)
addf.to_csv(out)
