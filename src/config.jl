@option struct SpecieConfig
    label::Symbol=:A
    "Can be: :Segment, :Small"
    type::Symbol=:Segment
    b::Float64=1.0
    M::Union{Float64, Nothing}=nothing
end

@option struct BlockConfig
    label::Symbol=:A
    segment::Symbol=:A
    length::Float64=0.5
    ends::Vector{Symbol}=[:AB]
end

@option struct ComponentConfig
    "Can be: :BCP, :SMOL"
    type::Symbol=:BCP
    label::Symbol=:AB
    length::Float64=1.0
    volume_fraction::Float64=1.0
    blocks::Vector{BlockConfig}=[BlockConfig(), BlockConfig(; label=:B, segment=:B)]
end

@option struct PolymerSystemConfig
    label::Union{String, Nothing}=nothing
    species::Vector{SpecieConfig}=[SpecieConfig(), SpecieConfig(; label=:B)]
    χN_map::Vector{Vector{Any}}=[[:A, :B, 20.0]]
    components::Vector{ComponentConfig}=[ComponentConfig()]
    chain_density::Float64=1.0
end

Configurations.from_dict(::Type{SpecieConfig}, ::Type{Symbol}, s) = Symbol(s)
Configurations.from_dict(::Type{ComponentConfig}, ::Type{Symbol}, s) = Symbol(s)
Configurations.from_dict(::Type{BlockConfig}, ::Type{Symbol}, s) = Symbol(s)

function load_config(yamlfile, T=PolymerSystemConfig)
    d = YAML.load_file(yamlfile; dicttype=Dict{String, Any})
    return from_dict(T, d)
end

function save_config(yamlfile, config)
    d = to_dict(config, YAMLStyle)
    return YAML.write_file(yamlfile, d)
end