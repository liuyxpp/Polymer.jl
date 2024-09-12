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

@option struct BlockCopolymerConfig
    label::Symbol=:BCP
    species::Vector{SpecieConfig}=SpecieConfig[]
    blocks::Vector{BlockConfig}=BlockConfig[]
end

@option struct ComponentConfig
    "Can be: :BCP, :SMOL"
    type::Symbol=:BCP
    label::Symbol=:AB
    length::Float64=1.0
    volume_fraction::Float64=1.0
    blocks::Vector{BlockConfig}=BlockConfig[]
end

@option struct PolymerSystemConfig
    label::Union{String, Nothing}=nothing
    species::Vector{SpecieConfig}=[SpecieConfig(), SpecieConfig(; label=:B)]
    χN_map::Vector{Vector{Any}}=[[:A, :B, 20.0]]
    components::Vector{ComponentConfig}=[ComponentConfig(; blocks=[BlockConfig(), BlockConfig(; label=:B, segment=:B)])]
    chain_density::Float64=1.0
end

Configurations.from_dict(::Type{SpecieConfig}, ::Type{Symbol}, s) = Symbol(s)
Configurations.from_dict(::Type{ComponentConfig}, ::Type{Symbol}, s) = Symbol(s)
Configurations.from_dict(::Type{BlockConfig}, ::Type{Symbol}, s) = Symbol(s)
Configurations.from_dict(::Type{BlockCopolymerConfig}, ::Type{Symbol}, s) = Symbol(s)

"""
    load_config(yamlfile, T=PolymerSystemConfig; top=nothing)

Load a config of type `T` from a YAML file.

# Arguments

- `yamlfile::String`: The path to the yaml file.
- `T::Type`: The type of the config to be loaded.

# Keyword Arguments

- `top` is the key of the top level of the config in the yaml. If `top` is `nothing`, the whole yaml file will be loaded. Otherwise, only the part under the key `top` will be loaded.
"""
function load_config(yamlfile, T=PolymerSystemConfig; top=nothing)
    d = YAML.load_file(yamlfile; dicttype=Dict{String, Any})
    d = isnothing(top) ? d : d[top]
    return from_dict(T, d)
end

"""
    save_config(yamlfile, config)

Save a config to a YAML file.

# Arguments

- `yamlfile::String`: The path to the yaml file.
- `config`: The config to be saved.
"""
function save_config(yamlfile, config)
    d = to_dict(config, YAMLStyle)
    return YAML.write_file(yamlfile, d)
end