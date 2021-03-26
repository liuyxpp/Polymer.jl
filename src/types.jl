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

# traits for the type of chain architecture
abstract type PolymerArchitecture end
# Non-cyclic architectures
abstract type NonCyclicArchitecture <: PolymerArchitecture end
struct LinearArchitecture <: NonCyclicArchitecture end
abstract type BranchedArchitecture <: NonCyclicArchitecture end
struct StarArchitecture <: BranchedArchitecture end
struct CombArchitecture <: BranchedArchitecture end
struct GeneralBranchedArchitecture <: BranchedArchitecture end
# Polymer that has ring(s) in it.
abstract type CyclicArchitecture <: PolymerArchitecture end
struct RingArchitecture <: CyclicArchitecture end

iscyclicchain(::CyclicArchitecture) = true
iscyclicchain(::PolymerArchitecture) = false
isnoncyclicchain(pa::PolymerArchitecture) = !iscyclicchain(pa)

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
isfreeblockend(::BlockEnd) = false
isfreeblockend(::FreeEnd) = true

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

nblocks(bc::BlockCopolymer) = length(bc.blocks)

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
    χN_map::Union{Dict{Set{Symbol},AbstractFloat},Nothing}
    C::Real # = \rho_0 R_g^3 / N, dimensionless chain density

    function PolymerSystem(components::Vector{T}, conf, χN_map, C) where {T<:AbstractComponent}
        @argcheck _isasystem(components)
        @argcheck _isasystem(components, χN_map)
        new(components, conf, χN_map, C)
    end
end
PolymerSystem(components::Vector{T}; χN_map=nothing, conf=BulkConfinement(), C=1.0) where {T<:AbstractComponent} = PolymerSystem(components, conf, χN_map, C)

"""
Check if the volume fraction of all compnents sums to 1.0.
"""
function _isasystem(components)
    return mapreduce(x->x.ϕ, +, components) ≈ 1.0 ? true : false
end

isconfined(s::PolymerSystem) = isconfined(s.confinement)

multicomponent(s::PolymerSystem) = length(s.components) == 1 ? false : true
ncomponents(s::PolymerSystem) = length(s.components)

specie(s::AbstractSpecie) = s.label
specie(m::SmallMolecule) = m.label
specie(b::PolymerBlock) = specie(b.segment)

species(c::BlockCopolymer) = [specie(b) for b in c.blocks] |> unique |> sort
nspecies(c::BlockCopolymer) = species(c) |> length
species(m::SmallMolecule) = [specie(m)]
nspecies(m::SmallMolecule) = 1
function _species(components)
    sp = Symbol[]
    for c in components
        union!(sp, species(c))
    end
    return sp
end

species(c::Component) = species(c.molecule)
nspecies(c::Component) = species(c) |> length
species(s::PolymerSystem) = _species(s.components) |> sort
nspecies(s::PolymerSystem) = species(s) |> length

function _isasystem(components, χN_map)
    sp = _species(components)
    n = length(sp)
    # For one-specie system, always true.
    if n == 1
        return true
    end
    # for n > 1, `χN_map` should not be nothing
    if isnothing(χN_map)
        return false
    end
    # The number of key-val pairs should be equal to (n-1)n/2 with n being the number of species.
    if length(χN_map) != (n-1)*n÷2
        return false
    end

    # the species in all compnents should be consistent with the species in the χN_map
    sp2 = []
    for key in keys(χN_map)
        union!(sp2, key)
    end
    return Set(sp) == Set(sp2)
end

function systemtype(s::PolymerSystem)
    if !multicomponent(s)
        return NeatPolymer()
    end
    for c in s.components
        if c.molecule isa SmallMolecule
            return PolymerSolution()
        end
    end
    return PolymerBlend()
end