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

# traits for multi-species system
abstract type SpecieNumberType end
struct SingleSpecieSystem <: SpecieNumberType end
struct TwoSpeciesSystem <: SpecieNumberType end
struct MultiSpeciesSystem <: SpecieNumberType end

# traits for the type of charge distribution
abstract type ChargedType end
struct Neutral <: ChargedType end
struct SmearedCharge <: ChargedType end
struct DiscreteCharge <: ChargedType end

abstract type AbstractSpecie end
struct KuhnSegment{T} <: AbstractSpecie
    label::Symbol
    b::T # length
    M::T # molecular weight in g/mol
    κs::Vector{T}  # affinities to surfaces, can be empty. See Fredrison's book (2006) 193-194 for details

    KuhnSegment(label::Symbol, b::T, M::T, κs::Vector{T}) where T = new{T}(label, b, M, κs)
end

function KuhnSegment(label; b=1.0, M=1.0, κs=Float64[])
    if isempty(κs)
        b, M = promote(b, M)
        T = typeof(b)
        return KuhnSegment(Symbol(label), b, M, T[])
    else
        T = promote_type(typeof(b), typeof(M), eltype(κs))
        return KuhnSegment(Symbol(label), T(b), T(M), Vector{T}(κs))
    end
end

==(sp1::KuhnSegment, sp2::KuhnSegment) = sp1.label == sp2.label

Base.show(io::IO, s::KuhnSegment) = print(io, "KuhnSegment $(s.label) with b=$(s.b), M=$(s.M), κs=$(s.κs)")

abstract type BlockEnd end

struct FreeEnd <: BlockEnd
    label::Symbol
end
FreeEnd(; label=:EF) = FreeEnd(label)

struct BranchPoint <: BlockEnd
    label::Symbol
end

Base.show(io::IO, fe::FreeEnd) = print(io, "FreeEnd $(fe.label)")
Base.show(io::IO, bp::BranchPoint) = print(io, "BranchPoint $(bp.label)")

abstract type AbstractBlock end
struct PolymerBlock{T} <: AbstractBlock
    label::Symbol
    segment::KuhnSegment
    f::T # = N_b / N, N is the total number of segments in a whole polymer chain
    E1::BlockEnd
    E2::BlockEnd
end

function PolymerBlock(; label=:A, specie=:A, f=1.0, E1=FreeEnd(), E2=FreeEnd())
    return PolymerBlock(Symbol(label), KuhnSegment(specie), f, E1, E2)
end

function PolymerBlock(segment::KuhnSegment; label=segment.label, f=1.0, E1=FreeEnd(), E2=FreeEnd())
    return PolymerBlock(Symbol(label), segment, f, E1, E2)
end

Base.show(io::IO, b::PolymerBlock) = print(io, "PolymerBlock $(b.label) with f=$(b.f) of specie $(b.segment.label)")

function Base.show(io::IO, ::MIME"text/plain", b::PolymerBlock)
    f = round(b.f, digits=4)
    println(io, "PolymerBlock $(b.label) with f=$f of $(b.segment)")
    println(io, "    Chain ends:")
    println(io, "      * ", b.E1)
    print(io, "      * ", b.E2)
end

"""
- Check whether each block has a unqiue label.
- Check whether the length of all blocks in a chain sum to 1.0.
"""
function _isachain(blocks)
    labels = map(x->x.label, blocks)
    length(unique(labels)) == length(labels) || return false

    mapreduce(x->x.f, +, blocks) ≈ 1.0 || return false

    return true
end

abstract type AbstractMolecule end

struct SmallMolecule{T} <: AbstractMolecule
    label::Symbol
    b::T # length
    M::T # molecular weight in g/mol

    SmallMolecule(label::Symbol, b::T, M::T) where T = new{T}(label, b, M)
end

const SpecieUnion = Union{KuhnSegment, SmallMolecule}

SmallMolecule(label; b=1.0, M=1.0) = SmallMolecule(Symbol(label), promote(b, M)...)

==(sp1::SmallMolecule, sp2::SmallMolecule) = sp1.label == sp2.label

Base.show(io::IO, smol::SmallMolecule) = print(io, "SmallMolecule $(smol.label) with b=$(smol.b)")

abstract type AbstractPolymer <: AbstractMolecule end

struct BlockCopolymer{T<:AbstractBlock} <: AbstractPolymer
    label::Symbol
    blocks::Vector{T}

    function BlockCopolymer(label::Symbol, blocks::Vector{T}; check=true) where {T<:PolymerBlock}
        check && @argcheck _isachain(blocks)
        new{T}(label, blocks)
    end
end

function Base.show(io::IO, bcp::BlockCopolymer)
    n = length(bcp.blocks)
    println(io, "BlockCopolymer $(bcp.label) with $n blocks:")
    for b in bcp.blocks
        println(io, "  * $b")
    end
end

function Base.show(io::IO, ::MIME"text/plain", bcp::BlockCopolymer)
    n = length(bcp.blocks)
    println(io, "BlockCopolymer $(bcp.label) with $n blocks:")
    for b in bcp.blocks
        print(io, "  * ")
        show(io, "text/plain", b)
        println(io)
    end
end

struct RandomCopolymer <: AbstractPolymer end
struct AlternatingCopolymer <: AbstractPolymer end

struct Particle <: AbstractMolecule end
struct GiantMolecule <: AbstractMolecule end

abstract type AbstractComponent end

struct Component{T<:AbstractMolecule, S<:Real} <: AbstractComponent
    molecule::T
    α::S # N / N_ref, N_ref is the total number of segments in a reference polymer chain
    ϕ::S # = n*N/V, number density of this component in the system.

    function Component(molecule::T, α::S, ϕ::S) where {T<:AbstractMolecule, S}
        new{T, S}(molecule, α, ϕ)
    end
end

## Do not define `Base.show` for `Component` object to avoid disable tree view functionality of Pluto.jl.

# Base.show(io::IO, c::Component) = print(io, "Component $(c.molecule.label) with ϕ=$(c.ϕ) and α=$(c.α) contains $(c.molecule)")

# function Base.show(io::IO, ::MIME"text/plain", c::Component)
#     println(io, "Polymer system component $(c.molecule.label):")
#     print(io, "* ", c.molecule)
#     println(io, "* ϕ = $(c.ϕ)")
#     println(io, "* α = $(c.α)")
# end

Component(molecule::T; α=1.0, ϕ=1.0) where {T<:AbstractMolecule} = Component(molecule, promote(α, ϕ)...)

abstract type AbstractSystem end

"""
    PolymerSystem{T} <: AbstractSystem

The key of `χN_map` should be a two-element `Set`. Each element is the unique symbol for a specie. For example, the `χN_map` of an AB diblock copolymer is a `Dict` with one entry `Set([:A, :B]) => χN`.

Note: For Edwards Model A (homopolymer + implicit solvent), one has to add a dummy solvent component in the `components` array to make the function `multicomponent` return correct result.
"""
struct PolymerSystem{T} <: AbstractSystem
    components::Vector{Component}
    confinement::ConfinementType
    χNmatrix::χNMatrix{T}
    C::T # = \rho_0 R_g^3 / N, dimensionless chain density

    function PolymerSystem(components::Vector{S}, conf, χNmatrix::χNMatrix{T}, C) where {S<:AbstractComponent, T}
        @argcheck _isasystem(components)
        @argcheck _isasystem(components, χNmatrix)
        new{T}(components, conf, χNmatrix, T(C))
    end
end

function PolymerSystem(components::Vector{T}, χNmatrix::χNMatrix;
                    conf=BulkConfinement(), C=1.0) where {T<:AbstractComponent}
    return PolymerSystem(components, conf, χNmatrix, C)
end

"""
    PolymerSystem(components, χNmap; conf=BulkConfinement(), C=1.0)

Constructor for `PolymerSystem`. `components` is a collection of `Component` objects. `χNmap` is a dictionary of interaction pairs. See [`χNMatrix`](@ref) for how to write a valid `χNmap`.
"""
function PolymerSystem(components, χNmap; conf=BulkConfinement(), C=1.0)
    return PolymerSystem(collect(components), χNMatrix(χNmap); conf=conf, C=C)
end

## Here, we define a `Base.print` instead of a `Base.show` to show `PolymerSystem`. If we define `Base.show` for `PolymerSystem`, then the tree view of this object will be disabled and fall back to this `Base.show`.
# This a compromise to show `PoymerSystem` object in their best representation in both REPL and Pluto.jl.
# To get a clean plain text display of `PolymerSystem`, use `print(system)` where `system` is a `PolymerSystem` object.
function Base.print(io::IO, s::PolymerSystem)
    n = length(s.components)
    name = ""
    for i in 1:n
        label = string(s.components[i].molecule.label)
        name = i < n ? name * label * " + " : name * label
    end
    println(io, "PolymerSystem ($name) contains $n components:")
    println(io)
    for c in s.components
        print(io, "Component $(c.molecule.label) with ϕ=$(c.ϕ) and α=$(c.α) contains ")
        println(io, c.molecule)
    end
    print(io, "with ", s.χNmatrix)
    print(io, "and confined by ", s.confinement)
end

"""
Check if the volume fraction of all compnents sums to 1.0.
"""
function _isasystem(components)
    return mapreduce(x->x.ϕ, +, components) ≈ 1.0 ? true : false
end

function _isasystem(components, χNmatrix::χNMatrix)
    return _species(components) == species(χNmatrix)
end