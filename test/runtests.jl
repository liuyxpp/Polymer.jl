using Polymer
using Test

@testset "utils.jl" begin
    @test infinity(1.0) == Inf
    @test infinity(1.0; negative=true) == -Inf
    @test infinity(1) == typemax(1)
    @test infinity(1; negative=true) == typemin(1)

    @test infinity(Float64) == Inf
    @test infinity(Int) == typemax(Int)

    @test isinfinity(100) == false
    @test isinfinity(Inf) == true
    @test isinfinity(-Inf) == true
    @test isinfinity(0.0) == false
end
