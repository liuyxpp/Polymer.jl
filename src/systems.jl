## Convenient functions for creating polymer chains.

branchpoints(n, prefix="EB") = [BranchPoint(Symbol(prefix*string(i))) for i in 1:n]
freeends(n, prefix="A") = [FreeEnd(Symbol(prefix*string(i))) for i in 1:n]

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

function linearABC(fA=0.3, fB=0.3)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    sC = KuhnSegment(:C)
    eb = branchpoints(2)
    fe = freeends(2)
    A = PolymerBlock(:A, sA, fA, eb[1], fe[1])
    C = PolymerBlock(:C, sC, 1-fA-fB, eb[2], fe[2])
    B = PolymerBlock(:B, sB, fB, eb[1], eb[2])
    return BlockCopolymer(:ABC, [A, B, C])
end

## Convenient functions for creating non-polymer components.

solvent(; label=:S) = SmallMolecule(label)

## Convenient functions for creating polymer systems.

"Default AB diblock copolymer system: two species are A and B, lengths of both segmeents are 1.0. However, χN and fA can be changed by Keyword argument. Default values are: χN=20.0 and fA=0.5."
function AB_system(; χN=20.0, fA=0.5)
    polymer = Component(diblock_chain(; fA=fA))
    return PolymerSystem([polymer], Dict(Set([:A, :B])=>χN))
end

function ABC_system(; χABN=40.0, χACN=40.0, χBCN=40.0, fA=0.3, fB=0.4)
    polymer = Component(linearABC(fA, fB))
    return PolymerSystem([polymer],
                         Dict(
                             Set([:A,:B])=>χABN,
                             Set([:A,:C])=>χACN,
                             Set([:B,:C])=>χBCN
                             )
                        )
end

"A/B homopolymers binary blend."
function A_B_system(; χN=20.0, ϕA=0.5, αA=1.0, αB=1.0)
    polymerA = Component(homopolymer_chain(label=:A), αA, ϕA)
    polymerB = Component(homopolymer_chain(label=:B), αB, 1-ϕA)
    return PolymerSystem([polymerA, polymerB], Dict([:A, :B]=>χN))
end

"AB diblock copolymers / A homopolymers blend."
function AB_A_system(; χN=20.0, ϕAB=0.5, fA=0.5, α=0.5)
    polymerAB = Component(diblock_chain(; fA=fA), 1.0, ϕAB)
    polymerA = Component(homopolymer_chain(; label=:hA, segment=KuhnSegment(:A)), α, 1-ϕAB)
    return PolymerSystem([polymerAB, polymerA],
                         Dict(Set([:A, :B])=>χN))
end

"AB diblock copolymers + solvent solution."
function AB_S_system(; χNAB=20.0, χNAS=100.0, χNBS=100.0, ϕAB=0.1, fA=0.5, α=0.01)
    polymer = Component(diblock_chain(; fA=fA), 1.0, ϕAB)
    sol = Component(solvent(), α, 1-ϕAB)
    return PolymerSystem([polymer, sol],
                         Dict(
                             Set([:A,:B])=>χNAB,
                             Set([:A,:S])=>χNAS,
                             Set([:B,:S])=>χNBS
                             )
                        )
end

"A homopolymer + B homopolymer + solvent solution."
function A_B_S_system(; χNAB=20.0, χNAS=100.0, χNBS=100.0, ϕA=0.1, ϕB=0.1, αA=1.0, αB=1.0, αS=0.01)
    polymerA = Component(homopolymer_chain(; label=:A, segment=KuhnSegment(:A)), αA, ϕA)
    polymerB = Component(homopolymer_chain(; label=:B, segment=KuhnSegment(:B)), αB, ϕB)
    sol = Component(solvent(), αS, one(ϕA)-ϕA-ϕB)
    return PolymerSystem([polymerA, polymerB, sol],
                         Dict(
                             Set([:A,:B])=>χNAB,
                             Set([:A,:S])=>χNAS,
                             Set([:B,:S])=>χNBS
                             )
                        )
end

"A homopolymer + B homopolymer + solvent1 + solvent2 solution."
function A_B_S1_S2_system(; χNAB=20.0, χNAS1=100.0, χNBS1=100.0, χNAS2=100.0, χNBS2=100.0, χNS1S2=100.0, ϕA=0.1, ϕB=0.1, ϕS1=0.4, αA=1.0, αB=1.0, αS1=0.01, αS2=0.01)
    polymerA = Component(homopolymer_chain(; label=:A, segment=KuhnSegment(:A)), αA, ϕA)
    polymerB = Component(homopolymer_chain(; label=:B, segment=KuhnSegment(:B)), αB, ϕB)
    S1 = Component(solvent(label=:S1), αS1, ϕS1)
    S2 = Component(solvent(label=:S2), αS2, one(ϕA)-ϕA-ϕB-ϕS1)
    return PolymerSystem([polymerA, polymerB, S1, S2],
                         Dict(
                             Set([:A,:B])=>χNAB,
                             Set([:A,:S1])=>χNAS1,
                             Set([:B,:S1])=>χNBS1,
                             Set([:A,:S2])=>χNAS2,
                             Set([:B,:S2])=>χNBS2,
                             Set([:S1,:S2])=>χNS1S2
                             )
                        )
end