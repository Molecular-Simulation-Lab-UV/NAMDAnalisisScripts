#=
# Atom type
#
# Used to interface with .pdb files
#
# By Mateo Barria-Urenda
=#

module AtomType

export Atom, atom2line, dist, getnormal

using LinearAlgebra  #  Cross product to get a plane normal
import Base.Vector  # Type conversion

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

end # module
