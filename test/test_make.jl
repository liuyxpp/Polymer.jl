@testset "make.jl: make" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config = load_config(configfile, top="system")
    abs = make(config)
    @test systemtype(abs) == PolymerSolution()
end