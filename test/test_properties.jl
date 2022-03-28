using Polymer

@testset "properties.jl" begin
    model = AB_system()
    @test component_number_type(model) == MonoSystem()
    @test Polymer.ϕ̄(model, :A) == 0.5
    @test Polymer.ϕ̄(model, :B) == 0.5

    model = AB_A_system()
    bcp = Polymer.molecule(:AB, model)
    @test component_number_type(model) == BinarySystem()
    @test Polymer.ϕ̄(model, :A) == 0.75
    @test Polymer.ϕ̄(model, :B) == 0.25

    @test component_label(1, model) == :AB
    @test component_labels(model) == [:AB, :hA]
    @test component_id(:AB, model) == 1
    @test component_ids(model) == [1, 2]

    @test Polymer.ϕ(1, model) == 0.5
    @test Polymer.ϕ(:AB, model) == 0.5
    @test Polymer.ϕs(model) == [0.5, 0.5]

    @test Polymer.α(2, model) == 0.5
    @test Polymer.α(:hA, model) == 0.5
    @test Polymer.αs(model) == [1.0, 0.5]

    @test Polymer.molecule(1, model) == model.components[1].molecule
    @test Polymer.molecule(:AB, model) == model.components[1].molecule
    @test molecule_labels(model) == [:AB, :hA]
    @test molecule_label(1, model) == :AB
    @test molecule_label(Polymer.molecule(:AB, model)) == :AB
    @test molecule_ids(model) == [1, 2]
    @test molecule_id(:AB, model) == 1
    @test molecule_id(Polymer.molecule(:AB, model), model) == 1

    blockA = Polymer.block(:A, bcp)
    blockB = Polymer.block(:B, bcp)
    @test Polymer.blocks(bcp) == bcp.blocks
    @test block_labels(bcp) == [:A, :B]
    @test block_label(blockA) == :A
    @test block_label(1, bcp) == :A
    @test block_ids(bcp) == [1, 2]
    @test block_id(:A, bcp) == 1
    @test block_id(blockA, bcp) == 1
    @test Polymer.block(1, bcp) == blockA
    @test Polymer.block(:A, bcp) == Polymer.block(1, bcp)
    @test block_lengths(bcp) == [0.5, 0.5]
    @test block_length(1, bcp) == 0.5
    @test block_length(:A, bcp) == 0.5
    @test block_species(bcp) == [:A, :B]
    @test block_specie(1, bcp) == :A
    @test block_specie(:A, bcp) == :A
    @test block_bs(bcp) == [1.0, 1.0]
    @test block_b(1, bcp) == 1.0
    @test block_b(:A, bcp) == 1.0
    @test Polymer.b(:A, model) == 1.0
    @test Polymer.b(:B, model) == 1.0
    @test Polymer.bs(model) == [1.0, 1.0]

    model = AB_S_system()
    @test component_number_type(model) == BinarySystem()
    @test Polymer.ϕ̄(model, :A) == 0.05
    @test Polymer.ϕ̄(model, :B) == 0.05
    @test Polymer.ϕ̄(model, :S) == 0.9

    s3 = A_B_S_system()
    @test component_number_type(s3) == TernarySystem()

    s4 = A_B_S1_S2_system()
    @test component_number_type(s4) == MultiComponentSystem()
end

nothing