import prody
import numpy as np
from os import remove, path

def dynamic_wrap(pdb: prody.AtomGroup,
                 frame: prody.Frame,
                 reference_element: str,
                 wrapping_elements: list[str]) -> None:
    
    ref = prody.calcCenter(pdb.select(reference_element))
    elements = [pdb.select(e) for e in wrapping_elements]
    cell = frame.getUnitcell()
    for element in elements:
        current_element = prody.calcCenter(element)
        dist = ref - current_element
        if np.any(dist > cell[:3]/2):
            target = current_element - np.where(dist > cell[:3]/2, -cell[:3], np.zeros_like(dist))
            prody.moveAtoms(element, to=target)
        elif np.any(dist < cell[:3]/2):
            target = current_element - np.where(-dist > cell[:3]/2, cell[:3], np.zeros_like(dist))
            prody.moveAtoms(element, to=target)


def file_wrap(pdbName: str,
              dcdName: str,
              refName: str,
              selections: list[str],
              wrap: bool=True,
              wrap_selection: str='protein',
              outFile: str='output.dcd') -> None:

    pdb = prody.parsePDB(pdbName)
    dcd = prody.DCDFile(dcdName)

    dcd.link(pdb)
    dcd.setCoords(pdb)
    dcd.setAtoms(pdb.select(refName))

    if path.exists(outFile): remove(outFile)

    dcd2 = prody.DCDFile(outFile, mode='w')

    for frame in dcd:
        ref = prody.calcCenter(pdb.select(selections[0]))
        cell = frame.getUnitcell()
        for i in selections[1:]:
            current_monomer = prody.calcCenter(pdb.select(i))
            dist = ref - current_monomer
            if np.any(dist > cell[:3]/2):
                target = current_monomer - np.where(dist > cell[:3]/2, -cell[:3], np.zeros_like(dist))
                prody.moveAtoms(pdb.select(i), to=target)
            elif np.any(dist < cell[:3]/2):
                target = current_monomer - np.where(-dist > cell[:3]/2, cell[:3], np.zeros_like(dist))
                prody.moveAtoms(pdb.select(i), to=target)

        frame.superpose()
        if wrap:
            prody.wrapAtoms(pdb, unitcell=cell[:3], center=prody.calcCenter(pdb.select(wrap_selection)))

        dcd2.write(pdb, unitcell=cell)
    dcd2.close()