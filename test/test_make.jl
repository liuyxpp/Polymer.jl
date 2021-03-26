@testset "make.jl: load_config" begin
    configfile = joinpath(pwd(), "test/ABS.yml")
    config = load_config(configfile)
    @test config["System"]["label"] == "AB/S"
end

@testset "make.jl: make" begin
    configfile = joinpath(pwd(), "test/ABS.yml")
    config = load_config(configfile)
    abs = make(config)
    @test systemtype(abs) == PolymerSolution()
    chainAB = abs.components[1].molecule
    # @test chaintype(chainAB) == LinearArchitecture()
end