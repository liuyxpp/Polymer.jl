@testset "utils.jl" begin
    @test unicodesymbol2string('ϕ') == "phi"
    @test unicodesymbol2string('β') == "beta"
    @test unicodesymbol2string('b') == "b"
end