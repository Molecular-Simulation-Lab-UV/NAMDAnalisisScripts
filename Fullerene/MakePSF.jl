#!/bin/julia
desc = """
  Takes a pdb file of a fullerene and creates a topology.

  by Mateo Barria-Urenda
"""

using ArgParse  # To manage inputs

# Atom Type to work with pdb filer
include("./AtomType.jl")
using .AtomType

# Argument parsing

s = ArgParseSettings(description=desc)
@add_arg_table! s begin
    "file"
    help = "Input file"
    required = true
    "outname"
    help = "Output files (both a pdb and a psf)"
    required = true
end

pargs = parse_args(s)

# Store atoms
atoms = []
header = ""
# Read from file
open(pargs["file"]) do io
    global header
    lines = readlines(io)
    header = lines[1]
    for line in lines
        if length(line) >= 4 && line[1:4] == "ATOM"
            push!(atoms, Atom(line))
        end
    end
end

# Number of atoms
natoms = length(atoms)

## Find bonds
# Array where distances between atoms that are not carbons are set to infinite
C_distances = [[((strip(atoms[i].element) == strip(atoms[j].element) == "C")
                ? dist(atoms[i], atoms[j])
                : Inf)
                for j in 1:natoms] for i in 1:natoms]

# Array where distances between atoms that are carbons to oxygens are set to infinite
CO_distances = [[((strip(atoms[i].element) in ["O"]  && strip(atoms[j].element) in ["C"])
                ? dist(atoms[i], atoms[j])
                : Inf)
                for j in 1:natoms] for i in 1:natoms]

# Store all bonds
bonds = []
# List all bonds
for i in 1:natoms
    global bonds
    element = strip(atoms[i].element)
    if element == "C"
        # For C we find the closest three other C atoms
        closest_dist = sort(C_distances[i])[2:4]  # 1 is itself, 2:4 are the closest three atoms
        # Get inds
        inds = [findfirst(x -> x == closest_dist[j], C_distances[i]) for j in 1:3]
        # Get pairs of inds, sorted
        pairs = [sort([i, j]) for j in inds]
        for p in pairs
            if !(p in bonds)  # since they are sorted, we can check for repeats
                push!(bonds, p)
            end
        end
    elseif element == "O"
        # For O, we find the closest Carbon atom
        closest_dist = sort(CO_distances[i])[1]
        j = findfirst(x -> x == closest_dist, CO_distances[i])
        pair = sort([i, j])
        push!(bonds, pair)
    elseif element == "H"
        # For H we assume it's bonded to the previous atom
        push!(bonds, [i-1, i])
    end
end

sort!(bonds)

# Store all angles
angles = []
# Find angles by comparing all pairs of bonds
for i in 1:length(bonds)
    for j in i+1:length(bonds)
        A = bonds[i]
        B = bonds[j]
        if A[1] in B
            ind = findfirst(x -> x == A[1], B) # where is it repeated
            angle = [A[2], A[1], B[ind!=1 ? 1 : 2]]
            push!(angles, angle)
        elseif A[2] in B
            ind = findfirst(x -> x == A[2], B) # where is it repeated
            angle = [A[1], A[2], B[ind!=1 ? 1 : 2]]
            push!(angles, angle)
        end
    end
end

sort!(angles)

# Store dihedrals
dihedrals = []
# Find dihedrals by comparing all pairs of angles
for i in 1:length(angles)
    for j in i+1:length(angles)
        A = angles[i]
        B = angles[j]
        # check if there's a pair in common
        # since they are ordered, we don't need to check 1 and 3
        if any([strip(atoms[i].element) != "C" for i in [A ; B]])
            # only do dihedrals for Carbon
            continue
        elseif (A[1] == B[2] &&  A[2] in [B[1], B[3]])
            ind = findfirst(x -> x == A[2], B) # where is it repeated, either 1 or 3
            dihedral = [A[3], A[2], A[1], B[ind!=1 ? 1 : 3]]
        elseif (A[3] == B[2] &&  A[2] in [B[1], B[3]])
            ind = findfirst(x -> x == A[2], B) # where is it repeated, either 1 or 3
            dihedral = [A[1], A[2], A[3], B[ind!=1 ? 1 : 3]]
        else
            continue
        end
        if dihedral[end] < dihedral[1]  # flip
            reverse!(dihedral)
        end
        push!(dihedrals, dihedral)
    end
end

sort!(dihedrals)

println(bonds)
println(angles)
println(dihedrals)
