@testset "graph.jl: load_config" begin
    config = load_config("./ABS.yml")
    @test config["System"]["label"] == "AB/S"
end

@testset "graph.jl: make" begin
    config = load_config("./ABS.yml")
    abs = make(config)
    @test systemtype(abs) == PolymerSolution()
    chainAB = abs.components[1].molecule
    @test chaintype(chainAB) == LinearArchitecture()
end