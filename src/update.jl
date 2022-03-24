function update!(system::PolymerSystem, ϕs::AbstractVector, ::PolymerParameter{:ϕ})
    (sum(ϕs) == 1.0) || error("ϕs should have sum 1.0!")

    nc = ncomponents(system)
    (nc == 1) && return system

    (nc == length(ϕs)) || error("Length of ϕs should be equal to the number of components of the polymer system!")

    for i in 1:nc
        c = system.components[i]
        system.components[i] = Component(c.molecule, c.α, ϕs[i])
    end

    return system
end

function update!(system::PolymerSystem, ϕ::Real, ::PolymerParameter{:ϕ})
    nc = ncomponents(system)
    (nc == 2) || error("Binary system expected!")

    return update!(system, [ϕ, one(ϕ)-ϕ], ϕParam)
end