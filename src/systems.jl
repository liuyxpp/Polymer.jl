"""
Convenient functions for creating a polymer chain.
"""
function homopolymer_chain(; label=:A, segment=KuhnSegment(label), α=1.0, ϕ=1.0)
    labelE1 = Symbol(label, :1)
    labelE2 = Symbol(label, :2)
    A = PolymerBlock(label, segment, 1.0, FreeEnd(labelE1), FreeEnd(labelE2))
    return PolymerComponent(label, [A]; α=α, ϕ=ϕ)
end

function diblock_chain(; labelA=:A, labelB=:B, segmentA=KuhnSegment(labelA), segmentB=KuhnSegment(labelB), fA=0.5, fB=1-fA, α=1.0, ϕ=1.0)
    labelAB = Symbol(labelA, labelB)
    eAB = BranchPoint(labelAB)
    A = PolymerBlock(labelA, segmentA, fA, FreeEnd(labelA), eAB)
    B = PolymerBlock(labelB, segmentB, fB, FreeEnd(labelB), eAB)
    return PolymerComponent(labelAB, [A,B]; α=α, ϕ=ϕ)
end