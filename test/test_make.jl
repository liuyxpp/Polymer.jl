@testset "make.jl: load_config" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config = load_config(configfile)
    @test config isa PolymerSystemConfig
end

@testset "make.jl: make" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config = load_config(configfile)
    abs = make(config)
    @test systemtype(abs) == PolymerSolution()
end