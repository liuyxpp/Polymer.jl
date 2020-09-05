
## Convenient functions for creating polymer chains.

function homopolymer_chain(; label=:A, segment=KuhnSegment(label))
    labelE1 = Symbol(label, :1)
    labelE2 = Symbol(label, :2)
    A = PolymerBlock(label, segment, 1.0, FreeEnd(labelE1), FreeEnd(labelE2))
    return BlockCopolymer(label, [A])
end

function diblock_chain(; labelA=:A, labelB=:B, segmentA=KuhnSegment(labelA), segmentB=KuhnSegment(labelB), fA=0.5)
    labelAB = Symbol(labelA, labelB)
    fB = 1.0 - fA
    eAB = BranchPoint(labelAB)
    A = PolymerBlock(labelA, segmentA, fA, FreeEnd(labelA), eAB)
    B = PolymerBlock(labelB, segmentB, fB, FreeEnd(labelB), eAB)
    return BlockCopolymer(labelAB, [A,B])
end

## Convenient functions for creating non-polymer components.

solvent(; label=:S) = SmallMolecule(label)

## Convenient functions for creating polymer systems.

"AB diblock copolymers."
function AB_system()
    polymer = Component(diblock_chain())
    return PolymerSystem([polymer]; χN_map=Dict(Set([:A, :B])=>20.0))
end

"AB diblock copolymers / A homopolymers blend."
function AB_A_system()
    polymerAB = Component(diblock_chain(), 1.0, 0.5)
    polymerA = Component(homopolymer_chain(; label=:hA, segment=KuhnSegment(:A)), 0.5, 0.5)
    return PolymerSystem([polymerAB, polymerA];
                         χN_map=Dict(Set([:A, :B])=>20.0))
end

"AB diblock copolymers + solvent solution."
function AB_S_system()
    polymer = Component(diblock_chain(), 1.0, 0.1)
    sol = Component(solvent(), 0.01, 0.9)
    return PolymerSystem([polymer, sol];
                         χN_map=Dict(
                             Set([:A,:B])=>20.0,
                             Set([:A,:S])=>100.0,
                             Set([:B,:S])=>100.0
                             )
                        )
end