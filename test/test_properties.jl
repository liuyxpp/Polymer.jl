using Polymer
using Polymer:
    label, components, component, blocks, segment, block_ends, molecule

@testset "properties.jl" begin
    model = AB_system()
    @test component_number_type(model) == MonoSystem()
    @test specie_number_type(model) == TwoSpeciesSystem()
    @test Polymer.ϕ̄(model, :A) == 0.5
    @test Polymer.ϕ̄(model, :B) == 0.5

    model = AB_A_system()
    bcp = Polymer.molecule(:AB, model)
    @test component_number_type(model) == BinarySystem()
    @test specie_number_type(model) == TwoSpeciesSystem()
    @test Polymer.ϕ̄(model, :A) == 0.75
    @test Polymer.ϕ̄(model, :B) == 0.25

    @test component_label(1, model) == :AB
    @test component_labels(model) == [:AB, :hA]
    @test component_id(:AB, model) == 1
    @test component_ids(model) == [1, 2]

    @test Polymer.χN(model, :A, :B) == 20.0

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

    ϕc = ϕControlParameter(1, ϕParam)
    @test Polymer.getparam(model, ϕc) == 0.5
    αc = αControlParameter(2, αParam)
    @test Polymer.getparam(model, αc) == 0.5
    χNc = χNControlParameter(:A, :B, χNParam)
    @test Polymer.getparam(model, χNc) == 20.0
    fc = fControlParameter(1, 1, fParam)
    @test Polymer.getparam(model, fc) == 0.5
    bc = bControlParameter(:A, bParam)
    @test Polymer.getparam(model, bc) == 1.0

    model = AB_S_system()
    @test component_number_type(model) == BinarySystem()
    @test Polymer.ϕ̄(model, :A) == 0.05
    @test Polymer.ϕ̄(model, :B) == 0.05
    @test Polymer.ϕ̄(model, :S) == 0.9

    s3 = A_B_S_system()
    @test component_number_type(s3) == TernarySystem()
    @test specie_number_type(s3) == MultiSpeciesSystem()

    s4 = A_B_S1_S2_system()
    @test component_number_type(s4) == MultiComponentSystem()
    @test specie_number_type(s4) == MultiSpeciesSystem()
end

@testset "properties.jl: label" begin
    model = AB_A_system()
    cs = components(model)
    @test label.(cs) == component_labels(model)
    c1 = first(cs)
    @test label(c1) == :AB
    @test label(c1) == component_label(c1)
    AB = molecule(c1)
    @test label(AB) == :AB
    b1 = first(blocks(AB))
    @test label(b1) == :A
    s1 = segment(b1)
    @test label(s1) == :A
    e1 = first(block_ends(b1))
    @test label(e1) == :A

    @test label.(blocks(AB)) == block_labels(AB)
    @test label(1, model) == :AB
    @test isnothing(label(3, model))

    @test label(1, AB) == :A
    @test isnothing(label(0, AB))
end

@testset "properties.jl: specie, species" begin
    sol = solvent()
    @test length(specie_objects(sol)) == 1
    @test species(sol) == [:S]

    chain = homopolymer_chain()
    @test length(specie_objects(chain)) == 1
    @test species(chain) == [:A]

    chain = diblock_chain()
    @test length(specie_objects(chain)) == 2
    @test species(chain) == [:A, :B]

    chain = linearABC()
    @test length(specie_objects(chain)) == 3
    @test species(chain) == [:A, :B, :C]

    model = AB_system()
    @test length(specie_objects(model)) == 2
    @test species(model) == [:A, :B]

    model = AB_A_B_system()
    @test length(specie_objects(model)) == 2
    @test species(model) == [:A, :B]

    model = AB_S_system()
    @test length(specie_objects(model)) == 3
    @test species(model) == [:A, :B, :S]
    cs = components(model)
    c1 = first(cs)
    c2 = last(cs)
    @test species(c1) == [:A, :B]
    @test species(c2) == [:S]
    AB = molecule(c1)
    S = molecule(c2)
    @test species(AB) == [:A, :B]
    @test species(S) == [:S]
    bA = first(blocks(AB))
    bB = last(blocks(AB))
    @test specie(bA) == :A
    @test specie(bB) == :B
    @test specie(S) == :S
    @test isempty(blocks(S))

    model = AB_A_system()
    @test length(specie_objects(model)) == 2
    @test species(model) == [:A, :B]
end

nothing