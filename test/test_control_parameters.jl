@testset "control_parameters.jl: ControlParameter" begin
    fϕ(ϕ) = [0.8-ϕ, 0.2]  # TernarySystem, with third component ϕ fixed at 0.2.
    ϕc = ϕControlParameter(1, ϕParam, fϕ)
    @test ϕc(0.3) == [0.3, 0.5, 0.2]

    abs = A_B_S_system()
    ϕc2 = ϕControlParameter(:B, abs; func=fϕ)
    @test ϕc2(0.3) == [0.5, 0.3, 0.2]

    αc = αControlParameter(1, αParam)
    @test αc(1.1) == 1.1
    αc2 = αControlParameter(:A, abs)  # control α value for A component

    χNc = χNControlParameter(:A, :B, χNParam)
    @test χNc(10.0) == 10.0
    χNc2 = χNControlParameter(:A, :B, abs)

    ff(f) = [0.8-f, 0.2]  # triblock copolymer, with third block f fixed at 0.2
    fc = fControlParameter(1, 1, fParam, ff)
    @test fc(0.3) == [0.3, 0.5, 0.2]

    abc = ABC_system()
    fc2 = fControlParameter(1, :ABC, abc; func=ff)
    @test fc2(0.3) == [0.3, 0.5, 0.2]

    bc = bControlParameter(:A, bParam)
    @test bc(1.1) == 1.1
end

nothing