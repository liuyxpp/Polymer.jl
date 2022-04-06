specie_object(m::SmallMolecule) = m
specie_object(b::PolymerBlock) = b.segment

specie_objects(bcp::BlockCopolymer) = specie_object.(bcp.blocks) |> unique
specie_objects(m::SmallMolecule) = [m]
specie_objects(c::Component) = specie_objects(c.molecule)
function specie_objects(s::PolymerSystem)
    species = Any[]
    for c in s.components
        append!(species, specie_objects(c))
    end
    return unique(species)
end

function to_config(sp::KuhnSegment)
    return SpecieConfig(label=sp.label, type=:Segment, b=sp.b, M=sp.M)
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