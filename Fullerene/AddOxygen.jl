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
        default = .96
    "--angle",  "-a"
        help = "Angle of C-O-H in degrees"
        arg_type = Float64
        # default = 108.
        default = 130.
end

pargs = parse_args(s)

# Need at least indexes or random
if (length(pargs["indexes"]) < 1 && pargs["random"] == nothing)
    throw(ArgumentError("Must pass either --indexes or --random"))
end

angle = pargs["angle"]
angle = angle / 180 * π  # Convert to radians

"""
Atom(id, name, rescha, resid, x, y, z, occupancy, beta, atom)

  Type to parse atoms from a pdb file.
"""
struct Atom
    # TODO include all .pdb columns
    id::Int
    name::String
    # resname::String
    rescha::String
    resid::Int
    x::Float64
    y::Float64
    z::Float64
    occupancy::Float64
    beta::Float64
    atom::String
end

"""
Atom(line)

  Alternate constructor from a text line in .pdb format.
"""
function Atom(line::String)::Atom
    # TODO Use indexes instead of split to comply with format
    _, id, name, rescha, resid, x, y, z, occupancy, beta, atom = split(line)
    id = parse(Int, id)
    resid = parse(Int, resid)
    x = parse(Float64, x)
    y = parse(Float64, y)
    z = parse(Float64, z)
    occupancy = parse(Float64, occupancy)
    beta = parse(Float64, beta)
    return Atom(id, name, rescha, resid, x, y, z, occupancy, beta, atom)
end

"""
atom2line(Atom)

  Convert atom type back to a pdb file line.
"""
function atom2line(A::Atom)
    x, y, z, occupancy, beta = round.([A.x, A.y, A.z, A.occupancy, A.beta], digits=3)

    line = "ATOM"
    line *= " "^(7 - length("$(A.id)")) * "$(A.id)"
    line *= " "^(3 - length("$(A.name)")) * "$(A.name)"
    line *= " "^(8 - length("$(A.rescha)")) * "$(A.rescha)"
    line *= " "^(4 - length("$(A.resid)")) * "$(A.resid)"
    line *= " "^(12 - length("$(x)")) * "$(x)"
    line *= " "^(8 - length("$(y)")) * "$(y)"
    line *= " "^(8 - length("$(z)")) * "$(z)"
    line *= " "^(6 - length("$(occupancy)")) * "$(occupancy)"
    line *= " "^(6 - length("$(beta)")) * "$(beta)"
    line *= " "^(12 - length("$(A.atom)")) * "$(A.atom)"
    return line * "\n"
end

"""
dist(A, B)

  Get distance between two atoms
"""
function dist(A::Atom, B::Atom)::Float64
    x = B.x - A.x
    y = B.y - A.y
    z = B.z - A.z
    return sqrt(x^2 + y^2 + z^2)
end

"""
Vector(Atom)

  Convert an Atom type to a Vector with its coordinates (and only that)
"""
function Vector(A::Atom)
    return Vector([A.x, A.y, A.z])
end

"""
getnormal(A, B, C)

  Get normal vector describing the plane that passes through A, B, and C (either 3D vectors or Atoms)
"""
function getnormal(A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64})
    return normalize(cross((B - A), (C - A)))
end

function getnormal(A::Atom, B::Atom, C::Atom)
    return getnormal(Vector(A), Vector(B), Vector(C))
end


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
if  length(pargs["indexes"]) > 0
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
distances = [[dist(atoms[i], atoms[j])  for j in 1:natoms] for i in 1:natoms]

# Store all bonds
bonds = []
# Store new atoms
oxygens = []
hydrogens = []
index = natoms  # Will be used to numper oxygens and hydrogens
# for all atoms check closest three
for i in modify
    global index
    index += 1 # First oxygen is numbered after the last carbon.
    closest_dist = sort(distances[i])[2:4]  # 1 is itself, 2:4 are the closest three atoms
    inds = [findfirst(x -> x == closest_dist[j], distances[i]) for j in 1:3]
    normal = getnormal(atoms[inds]...)
    dir = dot(Vector(atoms[i]), normal)  # Compare normal to our carbon
    if dir < 0
        normal *= -1  # Flip normal in these cases
    end
    nr = norm(normal)  # length of normal vector, should be 1
    # angles of normal vector
    ntheta = acos(normal[3]/nr)
    nphi = sign(normal[2]) * acos(normal[1] / sqrt(normal[1]^2 + normal[2]^2))
    # Oxygen: just add normal vector * bond length to the atom
    x, y, z = [atoms[i].x, atoms[i].y, atoms[i].z] + normal * pargs["bondCO"]
    Oxygen = Atom(index, "O", "X", 1, x, y, z, atoms[i].occupancy, atoms[i].beta, "O")
    # Store
    push!(oxygens, Oxygen)
    index += 1
    # Hydrogen: Set it with angle relative to the normal vector
    new_theta = ntheta + (π - angle)
    x = x + (pargs["bondOH"] * sin(new_theta) * cos(nphi))
    y = y + (pargs["bondOH"] * sin(new_theta) * sin(nphi))
    z = z + (pargs["bondOH"] * cos(new_theta))
    Hydrogen = Atom(index, "H", "X", 1, x, y, z, atoms[i].occupancy, atoms[i].beta, "H")
    # Store
    push!(hydrogens, Hydrogen)
end

# Writeout
open(pargs["outfile"], "w") do io
    global header
    write(io, header * "\n")
    for atom in atoms
        write(io, atom2line(atom))
    end
    for i in 1:length(oxygens)
        write(io, atom2line(oxygens[i]))
        write(io, atom2line(hydrogens[i]))
    end
    write(io, "END")
end
