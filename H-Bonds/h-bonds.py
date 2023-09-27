import numpy
from datetime import datetime
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import \
HydrogenBondAnalysis as HBA
import argparse
import multiprocessing as mp
import time


parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through a pore.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-p', '--procs', type = int, required = False, help = 'Number of cores for the analysis')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'
trajStep = None
if arg.procs:
    nCPUs = arg.procs
else:
    nCPUs = 1

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0]: # Any amount of .dcd files. Multiple lines, in order up-down, can be input. Only requisite is "dcd" should appear as keyword.
        dcdName.append(l[1])
    elif l[0].lower() == 'psf':
        psfName = l[1]

    # I don't need wrapping yet. Will implement in the future
    # elif l[0].lower() == 'ref': # Alignment selection for superposition of the .dcd frames onto the pdb conformation.
    #     if len(l[1:]) > 1:
    #         refName = ' '.join(l[1:])
    #     else:
    #         refName = l[1]
    elif l[0].lower() == 'sel1': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            sel1Name = ' '.join(l[1:])
        else:
            sel1Name = l[1]
    elif l[0].lower() == 'sel2': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            sel2Name = ' '.join(l[1:])
        else:
            sel2Name = l[1]
    elif l[0].lower() == 'group1': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            group1 = ' '.join(l[1:])
        else:
            group1 = l[1]
    elif l[0].lower() == 'group2': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            group2 = ' '.join(l[1:])
        else:
            group2 = l[1]
    elif l[0].lower() == 'seq': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            trajStep = int(l[1])
    elif l[0].lower() == 'out': # Path to the output file, name included.
        outName = l[1]


def CalcHBonds(queue, universe, sel1=None, sel2=None, groups=None, trajSlice=[None, None, None]):
    """
    This function is a wrapper for HydrogenBondAnalysis to allow for parallelization via the
    multiprocessing Process class.

    Parameters
    ----------
    universe : MDAnalysis universe object, with psf and dcd information loaded.
    sel1 : First selection for hydrogen bonds calculation. E.g: protein.
        Default is `None`.
    sel2 : Equivalent to sel1. E.g: water. Hydrogens and acceptors are guessed for both sel1 and sel2.
        Default is `None`.
    groups : List (size 2) of interactions to be considered for the hydrogen bond analysis. E.g: [sel1, sel2]
        Default is `None`.
    trajSlice : List of size 3 with the start, stop and step-size to go over the trajectory.
        Default is [`None`, `None`, `None`] to go over all frames.

    If everything other than universe is `None`, then all hydrogen bonds are calculated for all frames.
    This is a computation intensive and long process, so be careful.

    Returns
    -------
    hbonds : MDAnalysis' HydrogenBondsAnalysis object, with all methods

    """

    hbonds = HBA(universe=universe, between=groups, update_selections=False)
    # Modify the following to accept *args and avoid the consecutive if statements.
    sel_acceptors = hbonds.guess_acceptors(sel1)
    sel_hydrogens = hbonds.guess_hydrogens(sel2)
    hbonds.hydrogens_sel = f"{sel_hydrogens}"
    hbonds.acceptors_sel = f"{sel_acceptors}"
    
    hbonds.run(start=trajSlice[0], stop=trajSlice[1], step=trajSlice[2], verbose=True)

    queue.put(hbonds.results.hbonds)

################ ----------------- ################

# Security check. Don't try to use more than the max amount of CPUs in the system.
if nCPUs > mp.cpu_count():
    nCPUs = mp.cpu_count()

def main():
    universes = []
    processes = []
    q = mp.Queue()
    for i in range(nCPUs):
        u = mda.Universe(psfName, dcdName)
        universes.append(u)
        iterStep = len(u.trajectory)/nCPUs
        if nCPUs - i == 1:
            trajSeq = [int(i*iterStep), len(u.trajectory), trajStep]
        else:
            trajSeq = [int(i*iterStep), int((i+1)*iterStep), trajStep]
        process = mp.Process(target=CalcHBonds, args=(q, universes[i], sel1Name, sel2Name, [group1, group2], trajSeq))
        processes.append(process)
        process.start()
    while q.qsize() < nCPUs:
        time.sleep(0.05)
    results = q.get()
    while q.empty() is not True:
        results = numpy.vstack((results, q.get()))
        for p in processes:
            p.join()
            time.sleep(0.1)
    results = results.astype('str') # Handy for writing to a file

    outFile = open(outName, 'w+')
    outFile.write('#Frame, Donor_index, Hydrogen_index, Acceptor_index, DA_distance, DHA_angle\n')

    for vals in results:
        outFile.write('{0} \n'.format(' \t '.join(vals)))
        
    outFile.close()


t1 = datetime.now()
if __name__ == '__main__':
    main()
t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
