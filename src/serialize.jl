function to_config(sp::KuhnSegment)
    κs = isempty(sp.κs) ? nothing : sp.κs
    return SpecieConfig(label=sp.label, type=:Segment, b=sp.b, M=sp.M, κs=κs)
end

function to_config(sp::SmallMolecule)
    return SpecieConfig(label=sp.label, type=:Small, b=sp.b, M=sp.M)
end

function to_config(b::PolymerBlock)
    ends = Symbol[]
    (b.E1 isa BranchPoint) && push!(ends, b.E1.label)
    (b.E2 isa BranchPoint) && push!(ends, b.E2.label)
    return BlockConfig(label=b.label, segment=specie(b), length=b.f, ends=ends)
end

function to_config(bcp::BlockCopolymer)
    sps = to_config.(specie_objects(bcp))
    bs = to_config.(blocks(bcp))
    return BlockCopolymerConfig(label=label(bcp), species=sps, blocks=bs)
end

function to_config(c::Component)
    type = c.molecule isa BlockCopolymer ? :BCP : :SMOL
    blocks = type == :BCP ? to_config.(c.molecule.blocks) : BlockConfig[]
    return ComponentConfig(type=type, label=c.molecule.label, length=c.α, volume_fraction=c.ϕ, blocks=blocks)
end

function to_config(χNmatrix::χNMatrix)
    map = []
    for ((sp1, sp2), v) in χNmatrix.map
        push!(map, [sp1, sp2, v])
    end
    return map
end

function to_config(system::PolymerSystem)
    species = to_config.(specie_objects(system))
    χN_map = to_config(system.χNmatrix)
    components = to_config.(system.components)
    return PolymerSystemConfig(species=species, χN_map=χN_map, components=components, chain_density=system.C)
end

"`from_config` is an alias for `make`"
const from_config = make