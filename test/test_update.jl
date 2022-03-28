using Polymer
import Polymer: update!, molecule, block_lengths, block_bs

function starAB3(; fA=0.4, fB1=0.2, fB2=0.2, fB3=0.2)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    A = PolymerBlock(:A, sA, fA, FreeEnd(:A), eb0)
    B1 = PolymerBlock(:B1, sB, fB1, FreeEnd(:B1), eb0)
    B2 = PolymerBlock(:B2, sB, fB2, FreeEnd(:B2), eb0)
    B3 = PolymerBlock(:B3, sB, fB3, FreeEnd(:B3), eb0)
    return BlockCopolymer(:AB3, [A, B1, B2, B3])
end

function AB3_A_system(; χN=12.0, ϕAB3=0.8, αAB3=1.0, αA=0.35)
    polymerAB3 = Component(starAB3(), αAB3, ϕAB3)
    polymerA = Component(homopolymer_chain(; label=:A, segment=KuhnSegment(:A)), αA, 1-ϕAB3)
    return PolymerSystem([polymerAB3, polymerA],
                         Dict(Set([:A, :B])=>χN))
end

@testset "update.jl: χNParam" begin
    system = AB3_A_system()
    update!(system, :A, :B, 10.0, χNParam)
    @test Polymer.χN(system, :A, :B) == 10.0

    system = AB3_A_system()
    update!(system, 11.0, χNParam)
    @test Polymer.χN(system, :A, :B) == 11.0

    system = AB3_A_system()
    update!(system, Dict(Set([:A, :B])=>13.0), χNParam)
    @test Polymer.χN(system, :A, :B) == 13.0
end

@testset "update.jl: molecule" begin
    system = AB3_A_system()
    bcp = starAB3(; fA=0.1, fB1=0.3, fB2=0.3, fB3=0.3)
    update!(system, 1, bcp)
    @test block_lengths(molecule(1, system)) == [0.1, 0.3, 0.3, 0.3]

    system = AB3_A_system()
    bcp = starAB3(; fA=0.1, fB1=0.3, fB2=0.3, fB3=0.3)
    update!(system, :AB3, bcp)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.3, 0.3, 0.3]

    system = AB3_A_system()
    bcp = starAB3(; fA=0.1, fB1=0.3, fB2=0.3, fB3=0.3)
    update!(system, molecule(:AB3, system), bcp)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.3, 0.3, 0.3]
end

@testset "update.jl: ϕParam" begin
    system = AB3_A_system()
    update!(system, [2, 1], [0.6, 0.4], ϕParam)
    @test Polymer.ϕ(1, system) == 0.4
    @test Polymer.ϕ(2, system) == 0.6

    system = AB3_A_system()
    update!(system, [:A, :AB3], [0.6, 0.4], ϕParam)
    @test Polymer.ϕ(:AB3, system) == 0.4
    @test Polymer.ϕ(:A, system) == 0.6

    system = AB3_A_system()
    update!(system, [0.4, 0.6], ϕParam)
    @test Polymer.ϕ(:AB3, system) == 0.4
    @test Polymer.ϕ(:A, system) == 0.6

    system = AB3_A_system()
    update!(system, 0.4, ϕParam)
    @test Polymer.ϕ(:AB3, system) == 0.4
    @test Polymer.ϕ(:A, system) == 0.6
end

@testset "update.jl: αParam" begin
    system = AB3_A_system()
    update!(system, 2, 0.5, αParam)
    @test Polymer.α(1, system) == 1.0
    @test Polymer.α(2, system) == 0.5

    system = AB3_A_system()
    update!(system, :A, 0.5, αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5

    system = AB3_A_system()
    update!(system, [2, 1], [0.5, 1.0], αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5

    system = AB3_A_system()
    update!(system, [:A, :AB3], [0.5, 1.0], αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5

    system = AB3_A_system()
    update!(system, [2], [0.5], αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5

    system = AB3_A_system()
    update!(system, [:A], [0.5], αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5

    system = AB3_A_system()
    update!(system, [1.0, 0.5], αParam)
    @test Polymer.α(:AB3, system) == 1.0
    @test Polymer.α(:A, system) == 0.5
end

@testset "update.jl: fParam" begin
    bcp = starAB3()
    update!(bcp, [2, 3, 4, 1], [0.2, 0.3, 0.4, 0.1], fParam)
    @test block_lengths(bcp) == [0.1, 0.2, 0.3, 0.4]

    bcp = starAB3()
    update!(bcp, [:B1, :B2, :B3, :A], [0.2, 0.3, 0.4, 0.1], fParam)
    @test block_lengths(bcp) == [0.1, 0.2, 0.3, 0.4]

    bcp = starAB3()
    update!(bcp, [0.1, 0.2, 0.3, 0.4], fParam)
    @test block_lengths(bcp) == [0.1, 0.2, 0.3, 0.4]

    bcp = diblock_chain()
    update!(bcp, 0.4, fParam)
    @test block_lengths(bcp) == [0.4, 0.6]

    system = AB3_A_system()
    update!(system, 1, [2, 3, 4, 1], [0.2, 0.3, 0.4, 0.1], fParam)
    @test block_lengths(molecule(1, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB3_A_system()
    update!(system, :AB3, [2, 3, 4, 1], [0.2, 0.3, 0.4, 0.1], fParam)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB3_A_system()
    update!(system, molecule(:AB3, system), [2, 3, 4, 1], [0.2, 0.3, 0.4, 0.1], fParam)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB3_A_system()
    update!(system, 1, [0.1, 0.2, 0.3, 0.4], fParam)
    @test block_lengths(molecule(1, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB3_A_system()
    update!(system, :AB3, [0.1, 0.2, 0.3, 0.4], fParam)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB3_A_system()
    update!(system, molecule(:AB3, system), [0.1, 0.2, 0.3, 0.4], fParam)
    @test block_lengths(molecule(:AB3, system)) == [0.1, 0.2, 0.3, 0.4]

    system = AB_A_system()
    update!(system, 1, 0.4, fParam)
    @test block_lengths(molecule(1, system)) == [0.4, 0.6]

    system = AB_A_system()
    update!(system, :AB, 0.4, fParam)
    @test block_lengths(molecule(:AB, system)) == [0.4, 0.6]

    system = AB_A_system()
    update!(system, molecule(:AB, system), 0.4, fParam)
    @test block_lengths(molecule(:AB, system)) == [0.4, 0.6]
end

@testset "update.jl: bParam" begin
    bcp = starAB3()
    Polymer._update!(bcp, 3, 1.2, bParam)
    @test block_bs(bcp) == [1.0, 1.0, 1.2, 1.0]

    bcp = starAB3()
    Polymer._update!(bcp, :B2, 1.2, bParam)
    @test block_bs(bcp) == [1.0, 1.0, 1.2, 1.0]

    bcp = starAB3()
    Polymer._update!(bcp, [:B3, :B1, :B2], [1.3, 1.1, 1.2], bParam)
    @test block_bs(bcp) == [1.0, 1.1, 1.2, 1.3]

    bcp = starAB3()
    Polymer._update!(bcp, [1.0, 1.2, 1.2, 1.2], bParam)
    @test block_bs(bcp) == [1.0, 1.2, 1.2, 1.2]

    system = AB3_A_system()
    update!(system, :A, 1.2, bParam)
    @test block_bs(molecule(1,system)) == [1.2, 1.0, 1.0, 1.0]
    @test block_bs(molecule(2,system)) == [1.2]

    system = AB3_A_system()
    update!(system, :B, 1.2, bParam)
    @test block_bs(molecule(1,system)) == [1.0, 1.2, 1.2, 1.2]
    @test block_bs(molecule(2,system)) == [1.0]

    system = AB3_A_system()
    update!(system, [:A, :B], [1.1, 1.2], bParam)
    @test block_bs(molecule(1,system)) == [1.1, 1.2, 1.2, 1.2]
    @test block_bs(molecule(2,system)) == [1.1]
    @test Polymer.b(:A, system) == 1.1
    @test Polymer.b(:B, system) == 1.2
end

nothing