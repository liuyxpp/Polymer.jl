
## Convenient functions for creating polymer chains.

function homopolymer_chain(; label=:A, segment=KuhnSegment(label), α=1.0, ϕ=1.0)
    labelE1 = Symbol(label, :1)
    labelE2 = Symbol(label, :2)
    A = PolymerBlock(label, segment, 1.0, FreeEnd(labelE1), FreeEnd(labelE2))
    return PolymerComponent(label, [A]; α=α, ϕ=ϕ)
end

function diblock_chain(; labelA=:A, labelB=:B, segmentA=KuhnSegment(labelA), segmentB=KuhnSegment(labelB), fA=0.5, α=1.0, ϕ=1.0)
    labelAB = Symbol(labelA, labelB)
    fB = 1.0 - fA
    eAB = BranchPoint(labelAB)
    A = PolymerBlock(labelA, segmentA, fA, FreeEnd(labelA), eAB)
    B = PolymerBlock(labelB, segmentB, fB, FreeEnd(labelB), eAB)
    return PolymerComponent(labelAB, [A,B]; α=α, ϕ=ϕ)
end

## Convenient functions for creating non-polymer components.

solvent(; label=:S, α=0.01, ϕ=0.0) = SmallMoleculeComponent(label; α=α, ϕ=ϕ)

## Convenient functions for creating polymer systems.

"AB diblock copolymers."
AB_system() = PolymerSystem([diblock_chain()])

"AB diblock copolymers / A homopolymers blend."
AB_A_system() = PolymerSystem([diblock_chain(; ϕ=0.5), homopolymer_chain(; label=:hA, ϕ=0.5)])

"AB diblock copolymers + solvent solution."
AB_S_system() = PolymerSystem([diblock_chain(; ϕ=0.5), solvent(; ϕ=0.5)])