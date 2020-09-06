import Polymer

@testset "utils.jl" begin
    @test unicodesymbol2string('ϕ') == "phi"
    @test unicodesymbol2string('β') == "beta"
    @test unicodesymbol2string('b') == "b"

    d = Dict("a"=>1, "b"=>2, "c"=>3)
    @test reverse_dict(d) == Dict(1=>"a", 2=>"b", 3=>"c")
    d = Dict("c"=>1, "a"=>1, "b"=>2)
    @test reverse_dict(d) == Dict(1=>"a", 2=>"b")
    d = Dict("a"=>1, "b"=>2, "c"=>1)
    @test reverse_dict(d) == Dict(1=>"a", 2=>"b")

    @test Polymer._sort_tuple2((2,1)) == (1,2)
    @test Polymer._sort_tuple2((1,2)) == (1,2)
end