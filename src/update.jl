
################
# update χNParam
################

function update!(system::PolymerSystem, sp1::Symbol, sp2::Symbol, χN, ::χNType)
    system.χNmatrix[sp1, sp2] = χN
    return system
end

function update!(system::PolymerSystem, χN::Real, ::χNType)
    (ncomponents(system) == 2) || return system
    sp1, sp2 = species(system)
    return update!(system, sp1, sp2, χN, χNParam)
end

function update!(system::PolymerSystem, χNmap, ::χNType)
    system.χNmatrix .= χNMatrix(χNmap)
    return system
end

#################
# update molecule
#################

function update!(system::PolymerSystem, id_mol::Integer, mol::AbstractMolecule)
    c = system.components[id_mol]
    system.components[id_mol] = @set c.molecule = mol
    return system
end

"""
`mol_old` can be a label or an instance of AbstractMolecule
"""
function update!(system::PolymerSystem, mol_old, mol_new::AbstractMolecule)
    return update!(system, molecule_id(mol_old, system), mol_new)
end

#################
# update ϕParam
#################

function update!(system::PolymerSystem, ids::AbstractVector{<:Integer}, ϕs::AbstractVector, ::ϕType)
    (sum(ϕs) == 1.0) || error("ϕs should have sum 1.0!")

    nc = ncomponents(system)
    (nc == 1) && return system

    (nc == length(ϕs)) || error("Length of ϕs should be equal to the number of components of the polymer system!")

    (nc == length(ids)) || error("Length of ids (or labels) should be equal to the number of components of the polymer system!")

    for i in 1:nc
        id = ids[i]
        c = system.components[id]
        system.components[id] = @set c.ϕ = ϕs[i]
    end

    return system
end

function update!(system::PolymerSystem, labels::AbstractVector{<:Symbol}, ϕs::AbstractVector, ::ϕType)
    ids = [component_id(label, system) for label in labels]
    return update!(system, ids, ϕs, ϕParam)
end

"""
The sequence of the input ϕs should be consistent with `component_ids`, i.e. ϕs[1] corresponds to the ϕ for components[1], etc.
"""
function update!(system::PolymerSystem, ϕs::AbstractVector, ::ϕType)
    return update!(system, 1:ncomponents(system), ϕs, ϕParam)
end

"""
Binary system is assumed implicitly.
"""
function update!(system::PolymerSystem, ϕ::Real, ::ϕType)
    nc = ncomponents(system)
    (nc == 2) || error("Binary system expected!")

    return update!(system, [ϕ, one(ϕ)-ϕ], ϕParam)
end

#################
# update αParam
#################

function update!(system::PolymerSystem, id::Integer, α::Real, ::αType)
    c = system.components[id]
    system.components[id] = @set c.α = α
    return system
end

function update!(system::PolymerSystem, label::Symbol, α::Real, ::αType)
    return update!(system, component_id(label, system), α, αParam)
end

function update!(system::PolymerSystem, ids::AbstractVector{<:Integer}, αs::AbstractVector, ::αType)
    nc = ncomponents(system)
    (nc == 1) && return system

    (length(ids) == length(αs)) || error("Length of ids and αs should be equal!")

    for (id, α) in zip(ids, αs)
        update!(system, id, α, αParam)
    end

    return system
end

function update!(system::PolymerSystem, labels::AbstractVector{<:Symbol}, αs::AbstractVector, ::αType)
    ids = [component_id(label, system) for label in labels]
    return update!(system, ids, αs, αParam)
end

"""
The sequence of the input ϕs should be consistent with `component_ids`, i.e. ϕs[1] corresponds to the ϕ for components[1], etc.
"""
function update!(system::PolymerSystem, αs::AbstractVector, ::αType)
    nc = ncomponents(system)
    (nc == length(αs)) || error("Length of αs should be equal to the number of components!")

    return update!(system, 1:nc, αs, αParam)
end

#################
# update fParam
#################

function update!(bcp::BlockCopolymer, ids::AbstractVector{<:Integer}, fs::AbstractVector, ::fType)
    nb = nblocks(bcp)
    (sum(fs) == 1.0) || error("All f in one blockcopolymer must = 1.")
    (nb == length(fs)) || error("Length of fs should be equal to the number of blocks of the block polymer!")
    (nb == length(ids)) || error("Length of ids (or labels) should be equal to the number of blocks of the block polymer!")

    for i in 1:length(fs)
        id = ids[i]
        b = block(id, bcp)
        bcp.blocks[id] = @set b.f = fs[i]
    end
    return bcp
end

function update!(bcp::BlockCopolymer, labels::AbstractVector{<:Symbol}, fs::AbstractVector, ::fType)
    ids = [block_id(label, bcp) for label in labels]
    return update!(bcp, ids, fs, fParam)
end

"""
The sequence of `fs` should be consistent with `block_ids`.
"""
function update!(bcp::BlockCopolymer, fs::AbstractVector, ::fType)
    return update!(bcp, 1:nblocks(bcp), fs, fParam)
end

"""
Two-block blockcopolymer and the first block is assumed.
"""
function update!(bcp::BlockCopolymer, f::Real, ::fType)
    (nblocks(bcp) == 2) || return bcp
    return update!(bcp, [f, one(f)-f], fParam)
end

function update!(system::PolymerSystem, id_mol::Integer, ids_or_labels::AbstractVector, fs::AbstractVector, ::fType)
    mol = molecule(id_mol, system)
    update!(mol, ids_or_labels, fs, fParam)
    return update!(system, id_mol, mol)
end

"""
`bcp` can be the label or instance of a BlockCopolymer.
"""
function update!(system::PolymerSystem, bcp, ids_or_labels::AbstractVector, fs::AbstractVector, ::fType)
    return update!(system, molecule_id(bcp, system), ids_or_labels, fs, fParam)
end

function update!(system::PolymerSystem, id_mol::Integer, fs::AbstractVector, ::fType)
    return update!(system, id_mol, 1:length(fs), fs, fParam)
end

"""
`bcp` can be the label or instance of a BlockCopolymer.
"""
function update!(system::PolymerSystem, bcp, fs::AbstractVector, ::fType)
    return update!(system, molecule_id(bcp, system), fs, fParam)
end

"""
Two-block block copolymer is assumed.
"""
function update!(system::PolymerSystem, id_mol::Integer, f::Real, ::fType)
    return update!(system, id_mol, [f, one(f)-f], fParam)
end

"""
Two-block block copolymer is assumed.
`bcp` can be the label or instance of a BlockCopolymer.
"""
function update!(system::PolymerSystem, bcp, f::Real, ::fType)
    return update!(system, molecule_id(bcp, system), f, fParam)
end

#################
# update bParam
#################

function update!(bcp::BlockCopolymer, id_block::Integer, b::Real, ::bType)
    bk = bcp.blocks[id_block]
    bcp.blocks[id_block] = @set bk.segment.b = b
    return bcp
end

function update!(bcp::BlockCopolymer, label_block::Symbol, b::Real, ::bType)
    id_block = block_id(label_block, bcp)
    return update!(bcp, id_block, b, bParam)
end

function update!(bcp::BlockCopolymer, ids::AbstractVector{<:Integer}, bs::AbstractVector, ::bType)
    (length(ids) == length(bs)) || error("Length of ids and bs should be equal!")

    for (id, b) in zip(ids, bs)
        update!(bcp, id, b, bParam)
    end
    return bcp
end

function update!(bcp::BlockCopolymer, labels::AbstractVector{<:Symbol}, bs::AbstractVector, ::bType)
    ids = [block_id(label, bcp) for label in labels]
    return update!(bcp, ids, bs, bParam)
end

"""
The sequence of `bs` should be consistent with `block_ids`.
"""
function update!(bcp::BlockCopolymer, bs::AbstractVector, ::bType)
    nb = nblocks(bcp)
    (nb == length(bs)) || error("Length of bs must be equal to the number of blocks!")

    return update!(bcp, 1:nb, bs, bParam)
end

"""
`id_mol` is the id of a BlockCopolymer.
`id_or_label` and `b` can be either scalar or vector.
"""
function update!(system::PolymerSystem, id_mol::Integer, id_or_label, b, ::bType)
    mol = molecule(id_mol, system)
    update!(mol, id_or_label, b, bParam)
    return update!(system, id_mol, mol)
end

"""
`bcp` is the label or instance of a BlockCopolymer.
`id_or_label` and `b` can be either scalar or vector.
"""
function update!(system::PolymerSystem, bcp, id_or_label, b, ::bType)
    return update!(system, molecule_id(bcp, system), id_or_label, b, bParam)
end

"""
The sequence of `bs` should be consistent with `block_ids`.
"""
function update!(system::PolymerSystem, id_mol::Integer, bs::AbstractVector, ::bType)
    return update!(system, id_mol, 1:length(bs), bs, bParam)
end

"""
`bcp` can be the label or instance of a BlockCopolymer.
The sequence of `bs` should be consistent with `block_ids`.
"""
function update!(system::PolymerSystem, bcp, bs::AbstractVector, ::bType)
    return update!(system, molecule_id(bcp, system), bs, bParam)
end
