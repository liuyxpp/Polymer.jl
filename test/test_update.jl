using Polymer

@testset "update.jl" begin
    system = AB_A_system()
    Polymer.update!(system, 0.4, ϕParam)
    @test system.components[1].ϕ == 0.4
    @test system.components[2].ϕ == 0.6
    Polymer.update!(system, [0.3, 0.7], ϕParam)
    @test system.components[1].ϕ == 0.3
    @test system.components[2].ϕ == 0.7
end