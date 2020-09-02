# traits for space dimension
abstract type SpaceDimension end
struct D1 <: SpaceDimension end
struct D2 <: SpaceDimension end
struct D3 <: SpaceDimension end

# traits for the type of confinement
abstract type ConfinementType end
struct BulkConfinement <: ConfinementType end
struct BoxConfinement <: ConfinementType end
struct SlabConfinement <: ConfinementType end
struct DiskConfinement <: ConfinementType end
struct SphereConfinement <: ConfinementType end
struct CylinderConfinement <: ConfinementType end
isconfined(::BulkConfinement) = false
isconfined(::ConfinementType) = true

# traits for the type of polymer system
abstract type PolymerSystemType end
struct NeatPolymer <: PolymerSystemType end
struct PolymerBlend <: PolymerSystemType end
struct PolymerSolution <: PolymerSystemType end

# traits for the type of polymer
abstract type PolymerType end
struct Homopolymer <: PolymerType end
abstract type Copolymer <: PolymerType end
struct BlockCopolymer <: Copolymer end
struct RandomCopolymer <: Copolymer end
struct AlternatingCopolymer <: Copolymer end

# traits for the type of chain architecture
abstract type PolymerArchitecture end
struct LinearArchitecture <: PolymerArchitecture end
abstract type BranchedArchitecture <: PolymerArchitecture end
struct StarArchitecture <: BranchedArchitecture end
struct CombArchitecture <: BranchedArchitecture end
# Polymer that has ring(s) in it.
struct RingArchitecture <: PolymerArchitecture end

"""
    islinearchain(<:PolymerArchitecture)

Check if the chain is linear. Currently we use explicit traits to check. It is possible to check the architecture by examining PolymerComponent instance. But it is not implemented yet.
"""
islinearchain(::LinearArchitecture) = true
islinearchain(::PolymerArchitecture) = false

# traits for the type of charge distribution
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
    KuhnSegment(label; b=1.0, M=1.0) = new(label, b, M)
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

abstract type AbstractBlock end
struct PolymerBlock <: AbstractBlock
    label::Symbol
    segment::KuhnSegment
    f::Real # = N_b / N, N is the total number of segments in a whole polymer chain
    E1::BlockEnd
    E2::BlockEnd
end

"""
Check if the length of all blocks in a chain sum to 1.0.
"""
function _isachain(blocks)
    mapreduce(x->x.f, +, blocks) == 1.0 ? true : false
end

abstract type AbstractComponent end
struct PolymerComponent{T<:AbstractBlock} <: AbstractComponent
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

    function SmallMoleculeComponent(label; specie=SmallMolecule(label), α=0.01, ϕ=0.0)
        new(label, specie, α, ϕ)
    end
end
struct ParticleComponent <: AbstractComponent end
struct GiantMoleculeComponent <: AbstractComponent end

"""
Note: For Edwards Model A (homopolymer + implicit solvent), one has to add a dummy solvent component in the `components` array to make the function `multicomponent` return correct result.
"""
struct PolymerSystem{T<:AbstractComponent}
    components::Vector{T}
    confinement::ConfinementType
    C::Real # = \rho_0 R_g^3 / N, dimensionless chain density

    function PolymerSystem(components::Vector{T}; conf=BulkConfinement(), C=1.0) where {T<:AbstractComponent}
        @argcheck _isasystem(components)
        new{T}(components, conf, C)
    end
end

"""
Check if the volume fraction of all compnents sums to 1.0.
"""
function _isasystem(components)
    mapreduce(x->x.ϕ, +, components) == 1.0 ? true : false
end

multicomponent(s::PolymerSystem) = length(s.components) == 1 ? false : true
ncomponents(s::PolymerSystem) = length(s.components)

species(c::PolymerComponent) = [b.segment.label for b in c.blocks] |> unique
nspecies(c::PolymerComponent) = species(c) |> length
species(c::SmallMoleculeComponent) = [c.specie.label]
nspecies(c::SmallMoleculeComponent) = 1
function species(s::PolymerSystem)
    sp = []
    for c in s.components
        union!(sp, species(c))
    end
    return sp
end
nspecies(s::PolymerSystem) = species(s) |> length

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