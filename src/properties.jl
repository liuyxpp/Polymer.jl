isconfined(s::PolymerSystem) = isconfined(s.confinement)

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