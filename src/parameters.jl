abstract type AbstractParameter end

"""
    PolymerParameter{T} <: AbstractParameter

Type represents phyiscal parameters of a polymer system. Each type with a specific T (normally a symbol) defines a parameter. E.g.

    ```julia
    PolymerParameter{:ϕ} # define the ϕ parameter
    ```

Pre-defined parameter symbols:
* χ: Flory-Huggins interaction parameters.
* N: chain length in number of monomers or Kuhn segments.
* f: volume fraction of a block in a block copolymer.
* ϕ: volume fraction of a type of chain in a polymer system.
* Rg: radius of gyration of a polymer chain.
* C: chain density C = ρ₀Rg^3/N.
* b: length of a monomer or a Kuhn segment.
* α: chain length of a chain with respect to a reference chain.
* τ: ratio of length of two or more blocks.

Fields
* `description`: a description of the parameter.
* `variable_name`: a string representation of the parameter which can be used as Julia variable name.
* `ascii_label`: a label with pure ascii characters which shall be used for file names, directory names, etc.
* `plot_label`: a LaTeXStrings used as label for plotting.
* `value_type`: the number type of the parameter, e.g. integer, float, complex, etc.
* `allowed_min`: minimum value allowed. If negative infinity, using `typemin(Type)`.
* `allowed_max`: maximum value allowed. If positive infinity, using `typemax(Type)`.
* `allowed_values`: only used for discrete finite set. If the allowed values are infinite, use `allowed_min` and `allowed_max` instead.
* `disallowed_values`: an array of values which are invalid.

User can define their own parameter symbols as long as they also define following two field values. Other fields are deduced automatically from parametric type `T`.
* `description`
* `plot_label`
"""
struct PolymerParameter{T, R<:Real} <: AbstractParameter
    description::String
    variable_name::String
    ascii_label::String
    plot_label::LaTeXString
    value_type::Type{R}
    allowed_min::Union{R, typeof(Inf)}
    allowed_max::Union{R, typeof(Inf)}
    allowed_values::AbstractVector{R}
    disallowed_values::AbstractVector{R}
    function PolymerParameter{T}(desc, plot_label, value_type::Type{R}; allowed_min=-Inf, allowed_max=Inf, allowed_values=R[], disallowed_values=R[]) where {T, R<:Real}
        varname = String(T)
        ascii_label = ""
        for c in varname
            ascii_label *= unicodesymbol2string(c)
        end
        new{T, R}(desc, varname, ascii_label, plot_label, value_type, allowed_min, allowed_max, allowed_values, disallowed_values)
    end
end

"""
A list of pre-defined polymer parameters.
"""
const χParam = PolymerParameter{:χ}(
    "Flory-Huggins Interaction Parameter.",
    L"\chi", Float64)
const NParam = PolymerParameter{:N}(
    "Chain length in number of monomers or Kuhn segments.",
    L"N", Int; allowed_min=zero(Int))
const χNParam = PolymerParameter{:χN}(
    "Flory-Huggins Interaction Parameter times N.",
    L"\chi N", Float64)
const fParam = PolymerParameter{:f}(
    "Volume fraction of a block in a block copolymer.",
    L"f", Float64; allowed_min=0, allowed_max=1)
const ϕParam = PolymerParameter{:f}(
    "volume fraction of a type of chain in a polymer system.",
    L"\phi", Float64; allowed_min=0, allowed_max=1)
const RgParam = PolymerParameter{:Rg}(
    "Radius of gyration of a polymer chain.",
    L"R_g", Float64; allowed_min=0)
const CParam = PolymerParameter{:C}(
    "Chain density C = ρ₀Rg^3/N.",
    L"C", Float64; allowed_min=0)
const bParam = PolymerParameter{:b}(
    "Length of a monomer or a Kuhn segment.",
    L"b", Float64; allowed_min=0)
const αParam = PolymerParameter{:α}(
    "Chain length of a chain with respect to a reference chain.",
    L"\alpha", Float64; allowed_min=0)
const τParam = PolymerParameter{:τ}(
    "Ratio of length of two or more blocks.",
    L"\tau", Float64; allowed_min=0)

"""
Accessors for the `PolymerParameter` type.
"""
description(p::AbstractParameter) = p.description
as_variable_name(p::AbstractParameter) = p.variable_name
as_ascii_label(p::AbstractParameter) = p.ascii_label
as_plot_label(p::AbstractParameter) = p.plot_label