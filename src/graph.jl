using LightGraphs
using TikzGraphs
using TikzPictures
import Base: show

struct BlockCopolymerGraph
    graph::SimpleGraph
    block2edge::Dict{PolymerBlock, Tuple{Int, Int}}
    edge2block::Dict{Tuple{Int, Int}, PolymerBlock}
    joint2node::Dict{BranchPoint, Int}
    node2joint::Dict{Int, BranchPoint}

    function BlockCopolymerGraph(c::BlockCopolymer)
        g, d1, d2 = build_graph(c)
        rd1 = reverse_dict(d1)
        rd2 = reverse_dict(d2)
        new(g, d1, rd1, d2, rd2)
    end
end

species(bcg::BlockCopolymerGraph) = [specie(b) for b in keys(bcg.block2edge)] |> unique

"""
    build_graph(c::BlockCopolymer)

Build the graph representation of a BlockCopolymer instance, i.e. a block polymer chain. A LightGraphs undirected Graph instance as well as two Dicts are returned.

The first Dict instance is a collection of `PolymerBlock => (Int, Int)` pair where each integer number is the vertex index in the Graph instance. Note that the order of these two integers in the tuple is important. The first one corresponds to the block E1 end and the second one to the block E2 end, respectively. By denfinition, it is a one-to-one mapping. Therefore, we can safely reverse this Dict to get another Dict of `(Int, Int) => PolymerBlock`.

The Second Dict instance is a collection of `BranchPoint => Int` pair. Similarly, the integer is also the node index in the Graph instance. By denfinition, it is also a one-to-one mapping. We can safely reverse this Dict to get another Dict of `Int => BranchPoint`.

A block polymer chain consists of one or more polymer blocks, each of which consists of a number of identical monomers. Polymer blocks are connected by convalent bonds. Each polymer block has two ends. When the end is not connected to other blocks, we call it a free end. Otherwise, we call it a branch point or a joint. We map the architecture of a block copolymer chain to an undirected graph. Free ends and branch points of blocks are nodes of the graph, and blocks themselves are edges of the graph. Therefore, a block is uniquely determined by two node ids.

We use LightGraphs package to describe the graph. Thus the node id is just an integer. The first added node id is 1.
"""
function build_graph(c::BlockCopolymer)
	e1 = 0 # the node id for the first block end.
    e2 = 0 # the node id for the second block end.
	g = Graph()
	d1 = Dict{PolymerBlock, Tuple{Int, Int}}()
	d2 = Dict{BranchPoint, Int}()
    for b in c.blocks
        if isfreeblockend(b.E1)
            # Each free end is a distinct node in the graph.
			add_vertex!(g)
			e1 = nv(g)
        else
            # For branch point, we have to make sure if it is already added to the graph.
            # If added, we just use its node id.
            # If not, we have to create a new node corresponds to this branch point.
			if b.E1 ∈ keys(d2)
				e1 = d2[b.E1]
			else
				add_vertex!(g)
				e1 = nv(g)
				d2[b.E1] = e1
			end
        end
        # Do the same thing for the second block end.
		if isfreeblockend(b.E2)
			add_vertex!(g)
			e2 = nv(g)
		else
			if b.E2 ∈ keys(d2)
				e2 = d2[b.E2]
			else
				add_vertex!(g)
				e2 = nv(g)
				d2[b.E2] = e2
			end
        end
        # A block always corresponds to a new edge in the graph.
		add_edge!(g, e1, e2)
		d1[b] = (e1, e2)
	end
	return g, d1, d2
end

"""
    node_styles(g::SimpleGraph)

Make each node of a graph invisible, eliminating its gap between itself and connecting edge(s).
"""
node_styles(bcg::BlockCopolymerGraph) = Dict(vertices(bcg.graph) .=> "empty nodes, minimum size=0pt, inner sep=0pt")

"""
    edge_labels(g::SimpleGraph, d)

Assign the length of a block to the label of its corresponding edge in the graph representation of the block copolymer chain.


"""
function edge_labels(bcg::BlockCopolymerGraph)
	od = Dict()
	for (k, v) in bcg.edge2block
		od[_sort_tuple2(k)] = v.f
	end
	return od
end

"""
    SPECIECOLORS

An array of LaTeX color strings for denoting the type of a monomer. Here we use Tableau_10 color scheme (Ref: https://jiffyclub.github.io/palettable/tableau/#tableau_10). The colors are re-ordered with the following squence:
[blue, red, green, gray, orange, violet, brown, pink, yellow, light blue].

This form of specifing color is found [here at StackExchange](https://tex.stackexchange.com/questions/24434/how-to-specify-a-fill-color-in-rgb-format-in-a-node-in-tikzpicture).
"""
const SPECIECOLORS = ["draw={rgb,255:red,31;green,119;blue,180}",
                      "draw={rgb,255:red,214;green,39;blue,40}",
                      "draw={rgb,255:red,44;green,160;blue,44}",
                      "draw={rgb,255:red,127;green,127;blue,127}",
                      "draw={rgb,255:red,255;green,127;blue,14}",
                      "draw={rgb,255:red,148;green,103;blue,189}",
                      "draw={rgb,255:red,140;green,86;blue,75}",
                      "draw={rgb,255:red,227;green,119;blue,194}",
                      "draw={rgb,255:red,188;green,189;blue,34}",
                      "draw={rgb,255:red,23;green,190;blue,207}"
                      ]

"""
    edge_styles(bcg::BlockCopolymerGraph; colors=nothing, bends=nothing)

Return a Dict of `edge => style`.

## Keywords
* `colors`: an array of LaTeX color strings in the order of the symbol representation of each monomer (specie) type. If number of elements in `colors` less than `n`, we use first `n` colors in SPECIECOLORS with `n` being the number of species in the block copolymer.
* `bends`: an array of bending directions (bool values) defined in PGF/Tikz graph macro. Left bending is given by `true` and right bending is givne by `false` in the array. The length of this array should be at least the number of blocks of a block copolymer (`nb`). If the length is larger, only first `nb` values are used. If the number of elements in `bends` is length than `nb`, we use an alternating array starting at `false`.
"""
function edge_styles(bcg::BlockCopolymerGraph; colors=[], bends=[])
    sps = sort(species(bcg))
	n = length(sps)
    if length(colors) < n
        colors = SPECIECOLORS
    end
    cd = Dict(sps .=> colors[1:n])

    nb = ne(bcg.graph) # number of polymer blocks
    if length(bends) < nb
        bends = [i%2==0 for i in 1:nb]
    end
    bends = [bend ? "left" : "right" for bend in reverse(bends)]

	od = Dict()
	for (k, v) in bcg.edge2block
		bend = pop!(bends)
		color = cd[specie(v)]
		style = "ultra thick, $color, bend $bend"
		od[_sort_tuple2(k)] = style
	end
	return od
end

"""
Show `BlockCopolymerGraph` object as an image.
"""
function Base.show(io::IO, mime::MIME"image/svg+xml", g::BlockCopolymerGraph)
    tp = TikzGraphs.plot(g.graph,
                    node_styles=node_styles(g),
					edge_labels=edge_labels(g),
					edge_styles=edge_styles(g),
					layout=Layouts.SpringElectrical(),
                    options="scale=3")
    show(io, mime, tp)
end

"""
Show `BlockCopolymer` object as an image.
"""
Base.show(io::IO, mime::MIME"image/svg+xml", g::BlockCopolymer) = show(io, mime, BlockCopolymerGraph(g))

"""
    plot_graph(g::BlockCopolymerGraph)
    plot_graph(g::BlockCopolymer)

Similar to `show` method. But this method can configure extra options for plotting the block copolymer chain. It also return a `TikzPicture` object.
"""
function plot_graph(g::BlockCopolymerGraph; colors=[], bends=[])
    return TikzGraphs.plot(g.graph,
                    node_styles=node_styles(g),
					edge_labels=edge_labels(g),
					edge_styles=edge_styles(g; colors=colors, bends=bends),
					layout=Layouts.SpringElectrical(),
                    options="scale=3")
end

plot_graph(g::BlockCopolymer; colors=[], bends=[]) = plot_graph(BlockCopolymerGraph(g); colors=colors, bends=bends)

"""
    save_graph(f::TikzPictures.SaveType, g::BlockCopolymerGraph)
    save_graph(f::TikzPictures.SaveType, g::BlockCopolymer)

Save `BlockCopolymer` or `BlockCopolymerGraph` object into a file. The format can be PDF, SVG, and LaTeX. `f` should be one of:

* `TikzPictures.PDF("/path/to/save/graph")`: output /path/to/save/graph.pdf.
* `TikzPictures.SVG("/path/to/save/graph")`: output /path/to/save/graph.svg.
* `TikzPictures.TEX("/path/to/save/graph")`: output /path/to/save/graph.tex.
"""
save_graph(f::TikzPictures.SaveType, g::BlockCopolymerGraph; colors=[], bends=[]) = TikzPictures.save(f, plot_graph(g; colors=colors, bends=bends))

save_graph(f::TikzPictures.SaveType, g::BlockCopolymer; colors=[], bends=[]) = save_graph(f, BlockCopolymerGraph(g); colors=colors, bends=bends)

"""
    chaintype(g::BlockCopolymer)
    chaintype(g::BlockCopolymerGraph)

Return the trait of `PolymerArchitecture`.

We know our chain is always a connected graph. Therefore, it is easy to check whether it has cycles in it. For acyclic connected graph, it is merely a tree. And a tree has exactly (n - 1) edges where n is the number of nodes (vertices).
"""
function chaintype(g::BlockCopolymerGraph)
    # Chain has cycle(s)/ring(s) if and only if its number of edges is different than n - 1 where n is the number of nodes in the graph.
    if ne(g.graph) != nv(g.graph) - 1
        return RingArchitecture()
    end
    ds = degree(g.graph)
    # Linear chain cannot have branch point or can have branch point(s) with degree less than or equal to 2.
    if maximum(ds) <= 2
        return LinearArchitecture()
    end
    # Comb chain should only have branch points of degree 3 in the backbone.
    if sum(ds .== 3) + sum(ds .== 1) == length(ds) && sum(ds .== 3) > 1
        return CombArchitecture()
    end
    # Star chain should have a unique branch point with maximum degree (n>2).
    # Other branch points (if any) should have maximum degree less than n.
    if sum(ds .== maximum(ds)) > 1
        return GeneralBranchedArchitecture()
    else
        return StarArchitecture()
    end
end

chaintype(bcp::BlockCopolymer) = chaintype(BlockCopolymerGraph(bcp))

islinearchain(bcp::BlockCopolymer) = islinearchain(chaintype(bcp))