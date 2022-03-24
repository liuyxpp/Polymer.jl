using Polymer

@testset "properties.jl" begin
    model = AB_system()
    @test component_number_type(model) == MonoSystem()
    @test Polymer.ϕ̄(model, :A) == 0.5
    @test Polymer.ϕ̄(model, :B) == 0.5

    model = AB_A_system()
    @test component_number_type(model) == BinarySystem()
    @test Polymer.ϕ̄(model, :A) == 0.75
    @test Polymer.ϕ̄(model, :B) == 0.25

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