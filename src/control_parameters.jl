"""
AbstractControlParameter provides a unified interface to represent a control parameter for a polymer system. The control parameter can be varied to change the properties of a polymer system.

Concrete types include:
* ϕControlParameter
* αControlParameter
* fControlParameter
* χNControlParameter
* bControlParameter

The functor of each concrete type will return either a list or a scalar values of the parameters which fully detemines its setting.
"""
abstract type AbstractControlParameter end

struct ϕControlParameter{F<:Function} <: AbstractControlParameter
    id::Integer  # the unique id for a molecule in a polymer system
    param::ϕType
    func::F
end

"""
Default function is for binary system.
"""
function ϕControlParameter(id::Integer, param::ϕType=ϕParam,
                           func=(ϕ)->[one(ϕ)-ϕ])
    ϕs = func(0.1)
    @argcheck sum(ϕs) + 0.1 == 1.0 "func should return a vector with a sum equal to 1!"

    return ϕControlParameter(id, param, func)
end

function ϕControlParameter(label::Symbol, s::PolymerSystem;
                           param::ϕType=ϕParam, func=(ϕ)->[one(ϕ)-ϕ])
    id = molecule_id(label, s)
    @argcheck !isnothing(id) "$label is not a component label!"

    return ϕControlParameter(id, param, func)
end

"""
This functor produce a vector of ϕ whose order is assumed to be consistent with the polymer system it acts on.
"""
(c::ϕControlParameter)(ϕ) = insert!(c.func(ϕ), c.id, ϕ)

struct αControlParameter <: AbstractControlParameter
    id::Integer  # the unique id for a molecule in a polymer system
    param::αType
end

αControlParameter(id, param::αType=αParam) = αControlParameter(id, param)

function αControlParameter(label::Symbol, s::PolymerSystem;
                           param::αType=αParam)
    id = molecule_id(label, s)
    @argcheck !isnothing(id) "$label is not a component label!"

    return αControlParameter(id, param)
end

(::αControlParameter)(α) = α

"""
Set([sp1, sp2]) will determine which χN is the control parameter. For example, sp1=:A, sp2=:B, χ_AB N will be the control parameter.
"""
struct χNControlParameter <: AbstractControlParameter
    sp1::Symbol  # the unique name for the specie of a polymer system
    sp2::Symbol  # the unique name for the specie of a polymer system
    param::χNType
end

χNControlParameter(sp1, sp2, param::χNType=χNParam) = χNControlParameter(sp1, sp2, param)

function χNControlParameter(sp1::Symbol, sp2::Symbol, s::PolymerSystem;
                            param::χNType=χNParam)
    @argcheck sp1 ∈ species(s) "$sp1 is not a specie in the polymer system!"
    @argcheck sp2 ∈ species(s) "$sp2 is not a specie in the polymer system!"

    return χNControlParameter(sp1, sp2, param)
end

(::χNControlParameter)(χN) = χN

struct fControlParameter{F<:Function} <: AbstractControlParameter
    id_block::Integer  # the unique id for a block in a molecule
    id_mol::Integer  # the unique id for a molecule in a polymer system
    param::fType
    func::F
end

"""
Default function is for diblock copolymer.
"""
function fControlParameter(id_block, id_mol, param::fType=fParam,
                           func=(f)->[one(f)-f])
    fs = func(0.1)
    @argcheck sum(fs) + 0.1 == 1.0 "func should return a vector with a sum equal to 1!"

    return fControlParameter(id_block, id_mol, param, func)
end

function fControlParameter(id_block::Integer, label::Symbol, s::PolymerSystem;
                           param::fType=fParam, func=(f)->[one(f)-f])
    id_mol = molecule_id(label, s)
    @argcheck !isnothing(id_mol) "$label is not a component label!"
    @argcheck 0 < id_block <= nblocks(molecule(id_mol, s)) "Polymer block #$id_block is not in the polymer component!"

    return fControlParameter(id_block, id_mol, param, func)
end

function fControlParameter(label_block::Symbol, label_mol::Symbol,
                           s::PolymerSystem;
                           param::fType=fParam, func=(f)->[one(f)-f])
    id_mol = molecule_id(label_mol, s)
    @argcheck !isnothing(id_mol) "$label_mol is not a component label!"

    id_block = block_id(label_block, molecule(id_mol, s))
    @argcheck !isnothing(id_block) "Polymer block $label_block is not in the polymer component!"

    return fControlParameter(id_block, id_mol, param, func)
end

"""
This functor produce a vector of f whose order is assumed to be consistent with the molecule it acts on.
"""
(c::fControlParameter)(f) = insert!(c.func(f), c.id_block, f)

struct bControlParameter <: AbstractControlParameter
    sp::Symbol  # the unique name for the specie of a polymer system
    param::bType
end

bControlParameter(sp, param::bType=bParam) = bControlParameter(sp, param)

function bControlParameter(sp::Symbol, s::PolymerSystem; param::bType=bParam)
    @argcheck sp ∈ species(s) "$sp is not a specie in the polymer system!"

    return bControlParameter(sp, param)
end

(::bControlParameter)(b) = b