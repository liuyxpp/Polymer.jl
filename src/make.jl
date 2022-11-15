using Random: randstring

"""
    load_config(yamlfile)

Load configurations which directs creation of a polymer system and its components from a YAML file given at `yamlfile`.

The configurations are stored in a nested Dicts.
"""
# load_config(yamlfile) = YAML.load_file(yamlfile)

"""
    ObjType{name}

A type representing the type of an object defined in Polymer.jl.

A `ObjType` object can be passed to the first argument of [`make`](@ref) to dispatch to the correction version of this method to create an object designated by `name`:

* :System: PolymerSystem
* :Component: Component
* :χN_map: Dict{Set{Symbol},AbstractFloat}
* :Specie: AbstractSpecie (KuhnSegment) and SmallMolecule
* :SMOL: SmallMolecule
* :BCP: BlockCopolymer
* :Block: PolymerBlock
* :BlockEnd: BlockEnd (FreeEnd and BranchPoint)
"""
struct ObjType{name} end

function make(::ObjType{:BlockEnd}, config)
    if isempty(config)
        return FreeEnd(Symbol(randstring(3))), FreeEnd(Symbol(randstring(3)))
    end
    if length(config) == 1
        return FreeEnd(Symbol(randstring(3))), BranchPoint(Symbol(config[1]))
    end
    e1, e2 = Symbol.(config)
    return BranchPoint(e1), BranchPoint(e2)
end

function make(::ObjType{:Block}, config, sps)
    label = Symbol(config["label"])
    if haskey(config, "segment") && !isnothing(config["segment"])
        slabel = Symbol(config["segment"])
    else
        slabel= label
    end
    segment = nothing
    for sp in sps
        if sp.label == slabel
            segment = sp
            continue
        end
    end
    f = config["length"]
    E1, E2 = make(ObjType{:BlockEnd}(), config["ends"])
    return PolymerBlock(label, segment, f, E1, E2)
end

function make(::ObjType{:BCP}, config, sps)
    label = Symbol(config["label"])
    blocks = [make(ObjType{:Block}(), c, sps) for c in config["blocks"]]
    return BlockCopolymer(label, blocks)
end

function make(::ObjType{:SMOL}, config, sps)
    label = Symbol(config["label"])
    for sp in sps
        if sp.label == label
            return sp
        end
    end
end

function make(::ObjType{:Component}, config, sps)
    molecule = make(ObjType{Symbol(config["type"])}(), config, sps)
    kwargs = Dict()
    if haskey(config, "length")
        kwargs[:α] = config["length"]
    end
    if haskey(config, "volume_fraction")
        kwargs[:ϕ] = config["volume_fraction"]
    end
    return Component(molecule; kwargs...)
end

function make(::ObjType{:Specie}, config)
    label = Symbol(config["label"])
    kwargs = Dict()
    if haskey(config, "b") && !isnothing(config["b"])
        kwargs[:b] = config["b"]
    end
    if haskey(config, "M") && !isnothing(config["M"])
        kwargs[:M] = config["M"]
    end
    if haskey(config, "type")
        t = config["type"]
    else
        t = "Segment"
    end
    if t == "Segment"
        return KuhnSegment(label; kwargs...)
    else
        return SmallMolecule(label; kwargs...)
    end
end

function make(::ObjType{:χN_map}, config)
    sp1, sp2, χN = first(config)
    χN_map = Dict{Set{Symbol}, eltype(χN)}()
    for a in config
        s1, s2, χN = a
        χN_map[Set([Symbol(s1), Symbol(s2)])] = χN
    end
    return χN_map
end

"""
    make(::ObjType{:System}, config)

Make an `PolymerSystem` object.

Note: `confinement` is not implmented.
"""
function make(::ObjType{:System}, config)
    sps = [make(ObjType{:Specie}(), c) for c in config["species"]]
    components = [make(ObjType{:Component}(), c, sps) for c in config["components"]]
    kwargs = Dict()
    if haskey(config, "χN_map")
        χN_map = make(ObjType{:χN_map}(), config["χN_map"])
    else
        error("Must provide an interaction map.")
    end
    if haskey(config, "chain_density")
        kwargs[:C] = config["chain_density"]
    end
    return PolymerSystem(components, χN_map; kwargs...)
end

make(config) = make(ObjType{:System}(), config)

function make(config::SpecieConfig)
    label = config.label
    kwargs = isnothing(config.M) ? (; b=config.b) : (; b=config.b, M=config.M)
    (config.type == :Segment) && return KuhnSegment(label; kwargs...)
    (config.type == :Small) && return SmallMolecule(label; kwargs...)
    error("Unknown specie type!")
end

function make(config::BlockConfig, sps)
    i = findfirst(sp -> sp.label == config.segment, sps)
    isnothing(i) && error("No specie found for block!")

    label = config.label
    if isempty(config.ends)
        E1, E2 = FreeEnd(Symbol(randstring(3))), FreeEnd(Symbol(randstring(3)))
    elseif length(config.ends) == 1
        E1, E2 = FreeEnd(Symbol(randstring(3))), BranchPoint(config.ends[1])
    else
        E1, E2 = BranchPoint(config.ends[1]), BranchPoint(config.ends[2])
    end

    return PolymerBlock(label, sps[i], config.length, E1, E2)
end

function make(config::ComponentConfig, sps)
    (config.type ∈ [:BCP, :SMOL]) || error("Only BCP and SMOL is allowed for molecule type!")
    (config.type == :BCP && isempty(config.blocks)) && error("At least one block is expected for a block copolymer!")

    if config.type == :SMOL
        i = findfirst(sp -> sp.label == config.label, sps)
        isnothing(i) && error("No specie found for small molecule found!")
        molecule = sps[i]
    else
        blocks = make.(config.blocks, Ref(sps))
        molecule = BlockCopolymer(config.label, blocks)
    end

    return Component(molecule; α=config.length, ϕ=config.volume_fraction)
end

function _make_χNmap(config)
    χN_map = Dict{Set{Symbol}, Float64}()
    for a in config
        s1, s2, χN = a
        χN_map[Set([Symbol(s1), Symbol(s2)])] = χN
    end
    return χN_map
end

function make(config::PolymerSystemConfig)
    sps = make.(config.species)
    components = make.(config.components, Ref(sps))
    χN_map = _make_χNmap(config.χN_map)
    return PolymerSystem(components, χN_map; C=config.chain_density)
end

function KuhnSegment(config::SpecieConfig)
    sp = make(config)
    (sp isa SmallMolecule) && error("Configuration is for SmallMolecule!")
    return sp
end

function SmallMolecule(config::SpecieConfig)
    sp = make(config)
    (sp isa KuhnSegment) && error("Configuration is for KuhnSegment!")
    return sp
end

function BlockCopolymer(config::ComponentConfig, sps)
    mol = make(config, sps).molecule
    (mol isa SmallMolecule) && error("Configuration is for SmallMolecule!")
    return mol
end

function BlockCopolymer(config::BlockCopolymerConfig)
    sps = make.(config.species)
    blocks = make.(config.blocks, Ref(sps))
    return BlockCopolymer(config.label, blocks)
end

function SmallMolecule(config::ComponentConfig, sps)
    mol = make(config, sps).molecule
    (mol isa BlockCopolymer) && error("Configuration is for BlockCopolymer")
    return mol
end

PolymerBlock(config::BlockConfig, sps) = make(config, sps)
Component(config::ComponentConfig, sps) = make(config, sps)
PolymerSystem(config::PolymerSystemConfig) = make(config)