#!/bin/julia
desc = """
  Takes a pdb file of an unsubstituted fullerene and adds oxygens at the requested points
  Must pass either a list of indexes with -i or an amount of OH to generate with -r.
  If both are passed, -i takes priority.

  Examples:

  \$julia AddOxygen.jl C60.pdb C60_OH24.pdb -r 24

  Takes an existing pdb C60.pdb and adds oxygens at 24 random positions.

  \$julia AddOxygen.jl C60.pdb C60_OH4.pdb -i 1 10 12 15

  Takes an existing pdb C60.pdb and adds oxygens at 4 specified positions.

  by Mateo Barria-Urenda
"""

using LinearAlgebra  #  Cross product to get a plane normal
using ArgParse  # To manage inputs
using Random  # Random angle for C-O-H and random distribution of -OH
import Base.Vector  # Type conversion

# Atom Type to work with pdb filer
include("./AtomType.jl")
using .AtomType

# Argument parsing

s = ArgParseSettings(description=desc)
@add_arg_table! s begin
    "file"
    help = "Input file"
    required = true
    "outfile"
    help = "Output file"
    required = true
    "--indexes", "-i"
    help = "Indexes (from 1 to N) of the carbons that will be bonded to -OH groups."
    nargs = '+'
    arg_type = Int64
    "--random", "-r"
    help = "Modify a random number of carbons instead of a given list."
    arg_type = Int64
    "--bondCO", "-O"
    help = "Length of C-O bond (in A)"
    arg_type = Float64
    default = 1.411
    "--bondOH", "-H"
    help = "Length of O-H bond (in A)"
    arg_type = Float64
    default = 0.96
    "--angle", "-a"
    help = "Angle of C-O-H in degrees"
    arg_type = Float64
    # default = 108.
    default = 130.0
end

pargs = parse_args(s)

# Need at least indexes or random
if (length(pargs["indexes"]) < 1 && pargs["random"] === nothing)
    throw(ArgumentError("Must pass either --indexes or --random"))
end

angle = pargs["angle"]
angle = angle / 180 * π  # Convert to radians

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

# Check which to modify
if length(pargs["indexes"]) > 0
    modify = pargs["indexes"]
    if any(modify .> natoms)
        throw(ArgumentError("At least one index is bigger than the number of atoms! ($natoms)"))
    end
else
    if pargs["random"] > natoms
        throw(ArgumentError("Requested random replacements is larger than number of atoms! ($natoms)"))
    end
    modify = shuffle(1:natoms)[1:pargs["random"]]
end
# Matrix of all distances
distances = [[dist(atoms[i], atoms[j]) for j in 1:natoms] for i in 1:natoms]

# Store all bonds
bonds = []
# Store new atoms
oxygens = []
hydrogens = []
serial = natoms  # Will be used to numper oxygens and hydrogens
# for all atoms check closest three
for i in modify
    # Data to store in atom
    global serial
    serial += 1 # First oxygen is numbered after the last carbon.
    altLoc = atoms[i].altLoc
    resName = atoms[i].resName
    chainID = atoms[i].chainID
    resSeq = atoms[i].resSeq
    iCode = atoms[i].iCode
    occupancy = 0.0
    beta = 0.0
    charge = atoms[i].charge
    segment = atoms[i].segment
    # Define plane
    closest_dist = sort(distances[i])[2:4]  # 1 is itself, 2:4 are the closest three atoms
    inds = [findfirst(x -> x == closest_dist[j], distances[i]) for j in 1:3]
    normal = getnormal(atoms[inds]...)
    dir = dot(Vector(atoms[i]), normal)  # Compare normal to our carbon
    if dir < 0
        normal *= -1  # Flip normal in these cases
    end
    nr = norm(normal)  # length of normal vector, should be 1
    # angles of normal vector
    ntheta = acos(normal[3] / nr)
    nphi = sign(normal[2]) * acos(normal[1] / sqrt(normal[1]^2 + normal[2]^2))
    # Oxygen: just add normal vector * bond length to the atom
    x, y, z = [atoms[i].x, atoms[i].y, atoms[i].z] + normal * pargs["bondCO"]
    Oxygen = Atom("ATOM", serial, "O", altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, beta, segment, "O", charge)
    # Store
    push!(oxygens, Oxygen)
    serial += 1
    # Hydrogen: Set it with angle relative to the normal vector
    new_theta = ntheta + (π - angle)
    x = x + (pargs["bondOH"] * sin(new_theta) * cos(nphi))
    y = y + (pargs["bondOH"] * sin(new_theta) * sin(nphi))
    z = z + (pargs["bondOH"] * cos(new_theta))
    Hydrogen = Atom("ATOM", serial, "H", altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, beta, segment, "H", charge)
    # Store
    push!(hydrogens, Hydrogen)
end

# Writeout
open(pargs["outfile"], "w") do io
    global header
    write(io, header * "\n")
    for atom in atoms
        write(io, atom2pdb(atom))
    end
    for i in 1:length(oxygens)
        write(io, atom2pdb(oxygens[i]))
        write(io, atom2pdb(hydrogens[i]))
    end
    write(io, "END")
end
