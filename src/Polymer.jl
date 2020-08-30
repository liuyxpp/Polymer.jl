module Polymer

using ArgCheck
using REPL: symbol_latex
using LaTeXStrings

include("utils.jl")
export infinity, isinfinity, unicodesymbol2string

include("parameters.jl")
export PolymerSystemType, NeatPolymer, PolymerBlend, PolymerSolution
export PolymerType, Homopolymer, Copolymer, BlockCopolymer, RandomCopolymer
export PolymerArchitecture, LinearArchitecture, BranchedArchitecture, StarArchitecture, CombArchitecture, RingArchitecture
export SpaceDimension
export ConfinementType, BulkConfinement, BoxConfinement, SlabConfinement, DiskConfinement, SphereConfinement, CylinderConfinement

export AbstractParameter, PolymerParameter
export χParam, NParam, χNParam, fParam, RgParam, CParam, bParam, αParam, τParam
export description, as_variable_name, as_ascii_label, as_plot_label

export ChargedType, Neutral, SmearedCharge, DiscreteCharge
export BlockEnd, FreeEnd, BranchPoint, PolymerBlock
export AbstractSpecie, KuhnSegment, SmallMolecule
export AbstractComponent, PolymerComponent, SmallMoleculeComponent, ParticleComponent, GiantMoleculeComponent, PolymerSystem
export islinearchain, isconfined, ischarged, multicomponent, systemtype

abstract type PolymerSystemType end
struct NeatPolymer <: PolymerSystemType end
struct PolymerBlend <: PolymerSystemType end
struct PolymerSolution <: PolymerSystemType end

abstract type PolymerType end
struct Homopolymer <: PolymerType end
abstract type Copolymer <: PolymerType end
struct BlockCopolymer <: Copolymer end
struct RandomCopolymer <: Copolymer end

abstract type PolymerArchitecture end
struct LinearArchitecture <: PolymerArchitecture end
abstract type BranchedArchitecture <: PolymerArchitecture end
struct StarArchitecture <: BranchedArchitecture end
struct CombArchitecture <: BranchedArchitecture end
# Polymer that has ring(s) in it.
struct RingArchitecture <: PolymerArchitecture end
islinearchain(::LinearArchitecture) = true
islinearchain(::PolymerArchitecture) = false

abstract type SpaceDimension end
struct D1 <: SpaceDimension end
struct D2 <: SpaceDimension end
struct D3 <: SpaceDimension end

abstract type ConfinementType end
struct BulkConfinement <: ConfinementType end
struct BoxConfinement <: ConfinementType end
struct SlabConfinement <: ConfinementType end
struct DiskConfinement <: ConfinementType end
struct SphereConfinement <: ConfinementType end
struct CylinderConfinement <: ConfinementType end
isconfined(::BulkConfinement) = false
isconfined(::ConfinementType) = true

abstract type ChargedType end
struct Neutral <: ChargedType end
struct SmearedCharge <: ChargedType end
struct DiscreteCharge <: ChargedType end
ischarged(::Neutral) = false
ischarged(::ChargedType) = true

abstract type AbstractSpecie end
struct KuhnSegment <: AbstractSpecie
    label::Symbol
    b::Real # length
    M::Real # molecular weight in g/mol
    KuhnSegment(label, b=1.0, M=1.0) = new(label, b, M)
end
struct SmallMolecule <: AbstractSpecie
    label::Symbol
    b::Real # length
    M::Real # molecular weight in g/mol
    SmallMolecule(label; b=1.0, M=1.0) = new(label, b, M)
end

abstract type BlockEnd end
struct FreeEnd <: BlockEnd
    label::Symbol
end
struct BranchPoint <: BlockEnd
    label::Symbol
end

struct PolymerBlock
    label::Symbol
    segment::KuhnSegment
    f::Real # = N_b / N, N is the total number of segments in a whole polymer chain
    E1::BlockEnd
    E2::BlockEnd
end

function _isachain(blocks)
    mapreduce(x->x.f, +, blocks) == 1.0 ? true : false
end
abstract type AbstractComponent end
struct PolymerComponent{T<:PolymerBlock} <: AbstractComponent
    label::Symbol
    blocks::Vector{T}
    architecture::PolymerArchitecture
    charged::ChargedType
    α::Real # N / N_ref, N_ref is the total number of segments in a reference polymer chain
    ϕ::Real # = n*N/V, number density of this component in the system.

    function PolymerComponent(label, blocks::Vector{T}; arch=LinearArchitecture(), charged=Neutral(), α=1.0, ϕ=1.0) where {T<:PolymerBlock}
        @argcheck _isachain(blocks)
        new{T}(label, blocks, arch, charged, α, ϕ)
    end
end
struct SmallMoleculeComponent <: AbstractComponent
    label::Symbol
    specie::SmallMolecule
    α::Real # 1 / N_ref
    ϕ::Real # n/V
end
struct ParticleComponent <: AbstractComponent end
struct GiantMoleculeComponent <: AbstractComponent end

function _isasystem(components)
    mapreduce(x->x.ϕ, +, components) == 1.0 ? true : false
end
"""
Note: For Edwards Model A (homopolymer + implicit solvent), one has to add a dummy solvent component in the `components` array to make the function `multicomponent` return correct result.
"""
struct PolymerSystem{T<:AbstractComponent}
    components::Vector{T}
    confinement::ConfinementType
    dim::SpaceDimension
    C::Real # = \rho_0 R_g^3 / N, dimensionless chain density

    function PolymerSystem(components::Vector{T}; conf=BulkConfinement(), dim=D1(), C=1.0) where {T<:AbstractComponent}
        @argcheck _isasystem(components)
        new{T}(components, conf, dim, C)
    end
end

multicomponent(s::PolymerSystem) = length(s.components) == 1 ? false : true

function systemtype(s::PolymerSystem)
    if !multicomponent(s)
        return NeatPolymer()
    end
    for c in s.components
        if c isa SmallMoleculeComponent
            return PolymerSolution()
        end
    end
    return PolymerBlend()
end

end # module
