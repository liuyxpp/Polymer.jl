"""
    load_config(yamlfile)

Load configurations which directs creation of a polymer system and its components from a YAML file given at `yamlfile`.

The configurations are stored in a nested Dicts.
"""
load_config(yamlfile) = YAML.load_file(yamlfile)

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
        return FreeEnd(), FreeEnd()
    end
    if length(config) == 1
        return FreeEnd(), BranchPoint(Symbol(config[1]))
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
    blocks = [make(ObjType{:Block}(), c, sps) for c in config["chain"]]
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
    χN_map = Dict()
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
    sps = [make(ObjType{:Specie}(), c) for c in config["Species"]]
    components = [make(ObjType{:Component}(), c, sps) for c in config["Components"]]
    kwargs = Dict()
    if haskey(config, "χN_map")
        kwargs[:χN_map] = make(ObjType{:χN_map}(), config["χN_map"])
    end
    if haskey(config, "chain_density")
        kwargs[:C] = config["chain_density"]
    end
    return PolymerSystem(components; kwargs...)
end

make(config) = make(ObjType{:System}(), config["System"])