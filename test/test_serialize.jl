@testset "serialize.jl: to_config" begin
    config_aba = to_config(AB_A_system())
    @test config_aba isa PolymerSystemConfig
    aba = make(config_aba)
    @test systemtype(aba) == PolymerBlend()
end

nothing