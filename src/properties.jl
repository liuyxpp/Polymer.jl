isconfined(::BulkConfinement) = false
isconfined(::ConfinementType) = true
isconfined(s::PolymerSystem) = isconfined(s.confinement)

ischarged(::Neutral) = false
ischarged(::ChargedType) = true

χN(s::PolymerSystem, sp1::Symbol, sp2::Symbol) = s.χNmatrix[sp1, sp2]
χNmap(s::PolymerSystem) = s.χNmatrix.map

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
    sps = Symbol[]
    for c in components
        append!(sps, species(c))
    end
    return unique(sps) |> sort
end

species(c::Component) = species(c.molecule)
nspecies(c::Component) = species(c) |> length
species(s::PolymerSystem) = _species(s.components)
nspecies(s::PolymerSystem) = species(s) |> length

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

function component_number_type(s::PolymerSystem)
    nc = ncomponents(s)
    if nc == 1
        t = MonoSystem()
    elseif nc == 2
        t = BinarySystem()
    elseif nc == 3
        t = TernarySystem()
    else
        t = MultiComponentSystem()
    end
    return t
end

function specie_number_type(s::PolymerSystem)
    ns = nspecies(s)
    if ns == 1
        t = SingleSpecieSystem()
    elseif ns == 2
        t = TwoSpeciesSystem()
    else
        t = MultiSpeciesSystem()
    end
    return t
end

component_label(c::Component) = c.molecule.label
component_labels(s::PolymerSystem) = [component_label(c) for c in s.components]

component_ids(s::PolymerSystem) = collect(1:ncomponents(s))
component_id(label::Symbol, s::PolymerSystem) = findfirst(component_labels(s) .== label)
component_id(c::Component, s::PolymerSystem) = component_id(component_label(c), s)

function component_label(id::Integer, s::PolymerSystem)
    (id < 1 || id > ncomponents(s)) && return nothing
    return component_labels(s)[id]
end

ϕs(s::PolymerSystem) = [c.ϕ for c in s.components]
ϕ(id::Integer, s::PolymerSystem) = ϕs(s)[id]
ϕ(label::Symbol, s::PolymerSystem) = ϕ(component_id(label, s), s)

αs(s::PolymerSystem) = [c.α for c in s.components]
α(id::Integer, s::PolymerSystem) = αs(s)[id]
α(label::Symbol, s::PolymerSystem) = α(component_id(label, s), s)

molecules(s::PolymerSystem) = [c.molecule for c in s.components]
molecule(id::Integer, s::PolymerSystem) = molecules(s)[id]
molecule(label::Symbol, s::PolymerSystem) = molecule(component_id(label, s), s)
molecule_labels(s::PolymerSystem) = component_labels(s)
molecule_label(id::Integer, s::PolymerSystem) = component_label(id, s)
molecule_label(mol::AbstractMolecule) = mol.label
molecule_ids(s::PolymerSystem) = component_ids(s)
molecule_id(label::Symbol, s::PolymerSystem) = component_id(label, s)
molecule_id(mol::AbstractMolecule, s::PolymerSystem) = molecule_id(molecule_label(mol), s)

isfreeblockend(::BlockEnd) = false
isfreeblockend(::FreeEnd) = true

nblocks(sm::SmallMolecule) = 0
nblocks(bc::BlockCopolymer) = length(bc.blocks)
nblocks(c::Component) = nblocks(c.molecule)
nblocks(s::PolymerSystem) = sum(nblocks.(s.components))

blocks(bcp::BlockCopolymer) = bcp.blocks
block_labels(bcp::BlockCopolymer) = [b.label for b in bcp.blocks]
block_label(b::AbstractBlock) = b.label
block_ids(bcp::BlockCopolymer) = collect(1:nblocks(bcp))
block_id(label::Symbol, bcp::BlockCopolymer) = findfirst(block_labels(bcp) .== label)
block_label(id::Integer, bcp::BlockCopolymer) = block_labels(bcp)[id]
block_id(b::AbstractBlock, bcp::BlockCopolymer) = block_id(block_label(b), bcp)

block(id::Integer, bcp::BlockCopolymer) = bcp.blocks[id]
block(label::Symbol, bcp::BlockCopolymer) = block(block_id(label, bcp), bcp)

block_lengths(bcp::BlockCopolymer) = [b.f for b in bcp.blocks]
block_length(id::Integer, bcp::BlockCopolymer) = block_lengths(bcp)[id]
block_length(label::Symbol, bcp::BlockCopolymer) = block_length(block_id(label, bcp), bcp)

block_species(bcp::BlockCopolymer) = [specie(b) for b in bcp.blocks]
block_specie(id::Integer, bcp::BlockCopolymer) = block_species(bcp)[id]
block_specie(label::Symbol, bcp::BlockCopolymer) = block_specie(block_id(label, bcp), bcp)

block_bs(bcp::BlockCopolymer) = [b.segment.b for b in bcp.blocks]
block_b(id::Integer, bcp::BlockCopolymer) = block_bs(bcp)[id]
block_b(label::Symbol, bcp::BlockCopolymer) = block_b(block_id(label, bcp), bcp)

function b(sp::Symbol, system::PolymerSystem)
    (sp ∈ species(system)) || error("$sp is not existing!")

    for mol in molecules(system)
        if mol isa SmallMolecule
            (sp == specie(mol)) && return mol.b
        else
            id_block = findfirst(sp .== block_species(mol))
            return block_b(id_block, mol)
        end
    end
end

bs(system) = [b(sp, system) for sp in species(system)]

getparam(system::PolymerSystem, cp::ϕControlParameter) = ϕ(cp.id, system)
getparam(system::PolymerSystem, cp::αControlParameter) = α(cp.id, system)
getparam(bcp::BlockCopolymer, cp::fControlParameter) = block_length(cp.id_block, bcp)
getparam(system::PolymerSystem, cp::fControlParameter) = getparam(molecule(cp.id_mol, system), cp)
getparam(system::PolymerSystem, cp::χNControlParameter) = χN(system, cp.sp1, cp.sp2)
getparam(system::PolymerSystem, cp::bControlParameter) = b(cp.sp, system)

"""
    ϕ̄(c::Component{SmallMolecule}, sp::Symbol)
    ϕ̄(c::Component{<:BlockCopolymer}, sp::Symbol)
    ϕ̄(s::PolymerSystem, sp::Symbol)

The averaging density (volume fraction) of specie `sp` in a `PolymerSystem` or a `Component`.

This functionality can be used to compute the enthalpy energy in the Flory-Huggins theory.
"""
function ϕ̄(c::Component{SmallMolecule}, sp::Symbol)
    (sp ∈ species(c)) || return 0.0
    return c.ϕ
end

function ϕ̄(c::Component{<:BlockCopolymer}, sp::Symbol)
    (sp ∈ species(c)) || return 0.0
    ϕ = 0.0
    for b in c.molecule.blocks
        (specie(b) == sp) && (ϕ += c.ϕ * b.f)
    end
    return ϕ
end

function ϕ̄(s::PolymerSystem, sp::Symbol)
    (sp ∈ species(s)) || return 0.0
    return sum([ϕ̄(c, sp) for c in s.components])
end