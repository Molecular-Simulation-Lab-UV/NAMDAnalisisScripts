#=
# Atom type
#
# Used to interface with .pdb files
#
# By Mateo Barria-Urenda
=#

module AtomType

export Atom, atom2pdb, dist, getnormal

using LinearAlgebra  #  Cross product to get a plane normal
import Base.Vector  # Type conversion

"""
Atom(id, name, rescha, resid, x, y, z, occupancy, beta, atom)

  Type to parse atoms from a pdb file.
"""
struct Atom
    record::String
    serial::Int
    name::String
    altLoc::String
    resName::String
    chainID::String
    resSeq::Int
    iCode::String
    x::Float64
    y::Float64
    z::Float64
    occupancy::Float64
    beta::Float64
    element::String
    charge::String
end

"""
Atom(line)

  Alternate constructor from a text line in .pdb format.
"""
function Atom(line::String)::Atom
    # try to access line, return " " if out of bounds
    function get_range(start::Int, stop::Int)
        if start > length(line) || stop > length(line)
            return " "
        else
            return line[start:stop]
        end
    end
    # try to access line, return " " if out of bounds
    function get_range(ind::Int)
        if ind > length(line)
            return " "
        else
            return string(line[ind])
        end
    end
    record = get_range(1, 6)
    serial = parse(Int, get_range(7, 11))
    name = get_range(13, 16)
    altLoc = get_range(17)
    resName = get_range(18, 20)
    chainID = get_range(22)
    resSeq = parse(Int, get_range(23, 26))
    iCode = get_range(27)
    x = parse(Float64, get_range(31, 38))
    y = parse(Float64, get_range(39, 46))
    z = parse(Float64, get_range(47, 54))
    occupancy = parse(Float64, get_range(55, 60))
    beta = parse(Float64, get_range(61, 66))
    element = get_range(77, 78)
    charge = get_range(79, 80)

    return Atom(record,
        serial,
        name,
        altLoc,
        resName,
        chainID,
        resSeq,
        iCode,
        x,
        y,
        z,
        occupancy,
        beta,
        element,
        charge)
end

"""
format_length(s, l, leading=true)

  Formats l to a strigth of length l.
  If leading=true, use leading spaces; otherwire use trailing spaces.

  If length(s) > l, it's cut to l.
  This is donefrom the front for Ints, and from the back for Strings and Floats.
"""
function format_length(s::AbstractString, l::Int, leading=true)
    if length(s) > l
        s = s[1:l]
    end
    if leading
        return " "^(l - length(s)) * s
    else
        return s * " "^(l - length(s))
    end
end

function format_length(s::Int, l::Int, leading=true)
    s = "$(s)"
    if length(s) > l
        s = s[end-l:end]
    end
    return format_length(s, l, leading)
end

function format_length(s::Float64, l::Int, leading=true)
    s = "$(s)"
    return format_length(s, l, leading)
end

"""
format_name(name)

  Name formatting has it's own rules.
"""
function format_name(name::AbstractString)
    if length(name) == 1
        " " * name * "  "
    elseif length(name) == 2
        name * "  "
    else
        format_length(name, 4, false)
    end
end

"""
atom2pdb(Atom)

  Convert atom type back to a pdb file line.
"""
function atom2pdb(A::Atom)
    # nshould be handle by format_length, but rounding is a bit better than truncating so
    x, y, z = round.([A.x, A.y, A.z], digits=3)
    occupancy, beta = round.([A.occupancy, A.beta], digits=2)

    line = format_length(A.record, 6, false)
    line *= format_length(A.serial, 5)
    line *= " "
    line *= format_name(A.name)
    line *= format_length(A.altLoc, 1)
    line *= format_length(A.resName, 3)
    line *= " "
    line *= format_length(A.chainID, 1)
    line *= format_length(A.resSeq, 4)
    line *= format_length(A.iCode, 1)
    line *= " "^3
    line *= format_length(x, 8)
    line *= format_length(y, 8)
    line *= format_length(z, 8)
    line *= format_length(occupancy, 6)
    line *= format_length(beta, 6)
    line *= " "^10
    line *= format_length(A.element, 2)
    line *= format_length(A.charge, 2)
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
