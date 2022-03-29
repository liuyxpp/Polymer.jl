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

@testset "parameters.jl: AbstractControlParameter" begin
    fϕ(ϕ) = [0.8-ϕ, 0.2]  # TernarySystem, with third component ϕ fixed at 0.2.
    ϕc = ϕControlParameter(1, ϕParam, fϕ)
    @test ϕc(0.3) == [0.3, 0.5, 0.2]

    αc = αControlParameter(1, αParam)
    @test αc(1.1) == 1.1

    χNc = χNControlParameter(:A, :B, χNParam)
    @test χNc(10.0) == 10.0

    ff(f) = [0.8-f, 0.2]  # triblock copolymer, with third block f fixed at 0.2
    fc = fControlParameter(1, 1, fParam, ff)
    @test fc(0.3) == [0.3, 0.5, 0.2]

    bc = bControlParameter(:A, bParam)
    @test bc(1.1) == 1.1
end

nothing