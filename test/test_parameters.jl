using Polymer

@testset "parameters.jl: AbstractParameter" begin
    @test ϕParam.allowed_min == 0.0
    @test ϕParam.allowed_max == 1.0
    @test χParam.allowed_min == -Inf
    @test χParam.allowed_max == Inf
    @test NParam.allowed_min == 0
    @test NParam.allowed_max == Inf
    @test fParam.allowed_min == 0.0
    @test fParam.allowed_max == 1.0
end

nothing