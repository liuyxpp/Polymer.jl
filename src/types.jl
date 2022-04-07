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

# traits for the type of polymer system
abstract type PolymerSystemType end
struct NeatPolymer <: PolymerSystemType end
struct PolymerBlend <: PolymerSystemType end
struct PolymerSolution <: PolymerSystemType end

# traits for multicomponent system
abstract type ComponentNumberType end
struct MonoSystem <: ComponentNumberType end
struct BinarySystem <: ComponentNumberType end
struct TernarySystem <: ComponentNumberType end
struct MultiComponentSystem <: ComponentNumberType end

# traits for the type of charge distribution
abstract type ChargedType end
struct Neutral <: ChargedType end
struct SmearedCharge <: ChargedType end
struct DiscreteCharge <: ChargedType end

abstract type AbstractSpecie end
struct KuhnSegment <: AbstractSpecie
    label::Symbol
    b::Real # length
    M::Real # molecular weight in g/mol
end
KuhnSegment(label; b=1.0, M=1.0) = KuhnSegment(label, b, M)

abstract type BlockEnd end

struct FreeEnd <: BlockEnd
    label::Symbol
end
FreeEnd(; label=:EF) = FreeEnd(label)

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

function PolymerBlock(; label=:A, specie=:A, f=1.0, E1=FreeEnd(), E2=FreeEnd())
    return PolymerBlock(label, KuhnSegment(specie), f, E1, E2)
end

"""
Check if the length of all blocks in a chain sum to 1.0.
"""
function _isachain(blocks)
    mapreduce(x->x.f, +, blocks) ≈ 1.0 ? true : false
end

abstract type AbstractMolecule end

struct SmallMolecule <: AbstractMolecule
    label::Symbol
    b::Real # length
    M::Real # molecular weight in g/mol
end
SmallMolecule(label; b=1.0, M=1.0) = SmallMolecule(label, b, M)

abstract type AbstractPolymer <: AbstractMolecule end

struct BlockCopolymer{T<:AbstractBlock} <: AbstractPolymer
    label::Symbol
    blocks::Vector{T}

    function BlockCopolymer(label, blocks::Vector{T}) where {T<:PolymerBlock}
        @argcheck _isachain(blocks)
        new{T}(label, blocks)
    end
end

struct RandomCopolymer <: AbstractPolymer end
struct AlternatingCopolymer <: AbstractPolymer end

struct Particle <: AbstractMolecule end
struct GiantMolecule <: AbstractMolecule end

abstract type AbstractComponent end

struct Component{T<:AbstractMolecule} <: AbstractComponent
    molecule::T
    α::Real # N / N_ref, N_ref is the total number of segments in a reference polymer chain
    ϕ::Real # = n*N/V, number density of this component in the system.

    function Component(molecule::T, α, ϕ) where {T<:AbstractMolecule}
        new{T}(molecule, α, ϕ)
    end
end
Component(molecule::T; α=1.0, ϕ=1.0) where {T<:AbstractMolecule} = Component(molecule, α, ϕ)

abstract type AbstractSystem end

"""
The key of `χN_map` should be a two-element `Set`. Each element is the unique symbol for a specie. For example, the `χN_map` of an AB diblock copolymer is a `Dict` with one entry `Set([:A, :B]) => χN`.

Note: For Edwards Model A (homopolymer + implicit solvent), one has to add a dummy solvent component in the `components` array to make the function `multicomponent` return correct result.
"""
struct PolymerSystem <: AbstractSystem
    components::Vector{AbstractComponent}
    confinement::ConfinementType
    χNmatrix::χNMatrix
    C::Real # = \rho_0 R_g^3 / N, dimensionless chain density

    function PolymerSystem(components::Vector{T}, conf, χNmatrix::χNMatrix, C) where {T<:AbstractComponent}
        @argcheck _isasystem(components)
        @argcheck _isasystem(components, χNmatrix)
        new(components, conf, χNmatrix, C)
    end
end

PolymerSystem(components::Vector{T}, χNmatrix::χNMatrix; conf=BulkConfinement(), C=1.0) where {T<:AbstractComponent} = PolymerSystem(components, conf, χNmatrix, C)

PolymerSystem(components::Vector{T}, χNmap; conf=BulkConfinement(), C=1.0) where {T<:AbstractComponent} = PolymerSystem(components, χNMatrix(χNmap); conf=conf, C=C)

"""
Check if the volume fraction of all compnents sums to 1.0.
"""
function _isasystem(components)
    return mapreduce(x->x.ϕ, +, components) ≈ 1.0 ? true : false
end

function _isasystem(components, χNmatrix::χNMatrix)
    return _species(components) == species(χNmatrix)
end