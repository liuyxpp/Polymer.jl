using Polymer
using Test

@testset "parameters.jl" begin
    @test χParam.allowed_min == -Inf
    @test χParam.allowed_max == Inf
    @test NParam.allowed_min == 0
    @test NParam.allowed_max == Inf
    @test fParam.allowed_min == 0.0
    @test fParam.allowed_max == 1.0
end
