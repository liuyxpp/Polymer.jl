χNmapAB = Dict(Set((:A,:B))=>20.0)
χNmapABS = Dict(Set((:A,:B))=>20.0, Set((:A, :S))=>80.0, Set((:B, :S))=>120.0)
χNmapABS_bad = Dict(Set((:A,:B))=>20.0, Set((:A, :S))=>80.0)
χNmapABS_singular = Dict(Set((:A,:B))=>20.0, Set((:A, :S))=>0.0, Set((:B, :S))=>20.0)

@testset "chiN.jl: _species" begin
    @test Polymer._species(χNmapAB) == [:A, :B]
    @test Polymer._species(χNmapABS) == [:A, :B, :S]
end

@testset "chiN.jl: _check_consistency" begin
    @test Polymer._check_consistency(χNmapAB)
    @test Polymer._check_consistency(χNmapABS)
    @test !Polymer._check_consistency(χNmapABS_bad)
end

@testset "chiN.jl: _χNmap_to_matrix" begin
    mat = Polymer._χNmap_to_matrix(χNmapAB)
    @test mat == [0 20.0; 20.0 0]
    mat = Polymer._χNmap_to_matrix(χNmapABS)
    @test mat == [0 20.0 80.0; 20.0 0 120.0; 80.0 120.0 0]

    @test_throws ErrorException Polymer._χNmap_to_matrix(χNmapABS_singular)
end

@testset "chiN.jl: χNMatrix" begin
    m = χNMatrix(χNmapAB)
    @test m == [0 20.0; 20.0 0]
    m = χNMatrix(χNmapABS)
    @test m == [0 20.0 80.0; 20.0 0 120.0; 80.0 120.0 0]
    @test size(m) == (3,3)

    @test_throws ErrorException χNMatrix(χNmapABS_bad)
    @test_throws ErrorException χNMatrix(χNmapABS_singular)
end

@testset "chiN.jl: getindex" begin
    m = χNMatrix(χNmapABS)
    @test m[1] == 0
    @test m[7] == 80.0
    @test m[1, 2] == 20.0
    @test m[2, 1] == 20.0
    @test m[1, 3] == 80.0
    @test m[2, 3] == 120.0
    @test m[:A] == 0
    @test m[:B, :B] == 0
    @test m[:B, :A] == 20.0
    @test m[:A, :S] == 80.0
    @test m[:S, :B] == 120.0
end

@testset "chiN.jl: setindex" begin
    m = χNMatrix(χNmapABS)
    m[1] = 10.0
    @test m[1] == 0
    @test m[:A] == 0
    @test m[:A, :A] == 0

    m[7] = 88.0
    @test m[7] == 88.0
    @test m[1, 3] == 88.0
    @test m[3, 1] == 88.0

    m[:A, :A] = 99.0
    @test m[1] == 0
    @test m[1, 1] == 0
    @test m[:A] == 0
    @test m[:A, :A] == 0

    m[:S, :B] = 102.0
    @test m[:B, :S] == 102.0
    @test m[2, 3] == 102.0
    @test m[3, 2] == 102.0
end
