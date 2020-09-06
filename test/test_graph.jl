using LightGraphs

@testset "graph.jl: BlockCopolymerGraph" begin
    chainAB = diblock_chain()
    g = BlockCopolymerGraph(chainAB)
    @test nv(g.graph) == 3
    @test ne(g.graph) == 2
    @test g.block2edge == Dict(chainAB.blocks[1]=>(1,2), chainAB.blocks[2]=>(3,2))
    @test g.edge2block == Dict((1,2)=>chainAB.blocks[1], (3,2)=>chainAB.blocks[2])
    @test g.joint2node == Dict(chainAB.blocks[1].E2=>2)
    @test g.node2joint == Dict(2=>chainAB.blocks[1].E2)

    @test sort(species(g)) == [:A, :B]
end

@testset "graph.jl: Methods" begin
    chainAB = diblock_chain()
    g = BlockCopolymerGraph(chainAB)
    @test length(node_styles(g)) == nv(g.graph)
    @test length(edge_labels(g)) == ne(g.graph)
    @test length(edge_styles(g)) == ne(g.graph)

    @test chaintype(chainAB) == LinearArchitecture()
end