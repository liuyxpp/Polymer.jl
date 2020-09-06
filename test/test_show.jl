### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ b6251b9a-ee70-11ea-2d4a-e56c3bcee76d
using Polymer

# ╔═╡ 80b5df48-ee76-11ea-1ce1-47c07d90a110
using LightGraphs

# ╔═╡ 88653aa4-ee76-11ea-1d32-59a8a26d2715
using TikzGraphs

# ╔═╡ 27d3f3c4-ee85-11ea-31bd-8956d26c9c9b
using TikzPictures

# ╔═╡ c090e90e-ee72-11ea-1b34-7b5c1ff84d24
begin
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	sC = KuhnSegment(:C)
	sD = KuhnSegment(:D)
	eb1 = BranchPoint(:EB1)
	eb2 = BranchPoint(:EB2)
	eb3 = BranchPoint(:EB3)
	A1 = PolymerBlock(:A1, sA, 0.2, FreeEnd(:A1), eb1)
	B = PolymerBlock(:B, sB, 0.3, eb1, eb2)
	A2 = PolymerBlock(:A2, sA, 0.3, eb2, eb3)
	C1 = PolymerBlock(:C1, sC, 0.1, eb3, FreeEnd(:C1))
	C2 = PolymerBlock(:C2, sC, 0.05, eb2, FreeEnd(:C2))
	D = PolymerBlock(:D, sD, 0.05, eb3, FreeEnd(:D))
	ab = BlockCopolymer(:ABCAC, [A1, B, A2, C1, C2, D]) 
end

# ╔═╡ de11d6ae-ee73-11ea-380f-451eea158e9c
n = nspecies(ab)

# ╔═╡ 0163aa4c-ee74-11ea-2f1f-b50463c30d24
Polymer.isfreeblockend(ab.blocks[1].E1)

# ╔═╡ 8707d6d8-ee73-11ea-1ce8-71fa50acb0bf
function nblockends(c)
	nf = 0
	bps = []
	for b in c.blocks
		if isfreeblockend(b.E1)
			nf += 1
		else
			push!(bps, b.E1)
		end
		if isfreeblockend(b.E2)
			nf += 1
		else
			push!(bps, b.E2)
		end
	end
	nbp = Set(bps) |> length
	return nf + nbp
end

# ╔═╡ df60ce48-f01c-11ea-0ccd-83995099b387
nblockends(ab) == nv(BlockCopolymerGraph(ab).graph)

# ╔═╡ 821487a0-ee7a-11ea-2f7f-1de469ce28df
function build_graph(c)
	e1 = 0
	e2 = 0
	g = Graph()
	d1 = Dict{PolymerBlock, Tuple{Int, Int}}()
	d2 = Dict{BranchPoint, Int}()
	for b in c.blocks
		if isfreeblockend(b.E1)
			add_vertex!(g)
			e1 = nv(g)
		else
			if b.E1 ∈ keys(d2)
				e1 = d2[b.E1]
			else
				add_vertex!(g)
				e1 = nv(g)
				d2[b.E1] = e1
			end
		end
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
		add_edge!(g, e1, e2)
		d1[b] = (e1, e2)
	end
	return g, d1, d2
end

# ╔═╡ 0bde3d3a-ee7f-11ea-18cc-f3ee75c77ce1
g, d1, d2 = build_graph(ab)

# ╔═╡ 06ee1f2c-ee83-11ea-3f3a-adec2515cc8f
rd1 = Dict(v => k for (k, v) in d1)

# ╔═╡ 2784cf78-ee83-11ea-363d-63a4fc2d0396
rd2 = Dict(v => k for (k, v) in d2)

# ╔═╡ 39653f5a-ee83-11ea-0560-4d7ab9eaf280
node_styles(g) = Dict(vertices(g) .=> "empty nodes, minimum size=0pt, inner sep=0pt")

# ╔═╡ eebf4d96-ee83-11ea-37ed-01cb453eec53
node_styles(g)

# ╔═╡ 67e66760-ee91-11ea-14cc-07e697ba072c
function _sort_tuple(t)
	a, b = t
	return a > b ? (b, a) : (a, b)
end

# ╔═╡ 9fab2686-ee91-11ea-106d-b3c9ca2e4231
_sort_tuple((3,2))

# ╔═╡ 040ebc18-ee84-11ea-163c-1fd9f5536e09
function edge_labels(g, d)
	od = Dict()
	for (k, v) in d
		od[_sort_tuple(k)] = v.f
	end
	return od
end

# ╔═╡ a8bc6e66-ee84-11ea-0fcf-71f636966091
edge_labels(g, rd1)

# ╔═╡ 07aaf17e-ee85-11ea-0723-677e9f5ae9e5
function edge_styles(g, d, c)
	bends = ["right", "left", "right", "left", "right", "left"]
	colors = ["draw={rgb,255:red,31;green,119;blue,180}", "draw={rgb,255:red,214;green,39;blue,40}", "draw={rgb,255:red,44;green,160;blue,44}", "draw={rgb,255:red,127;green,127;blue,127}"]
	sps = sort(species(c))
	n = length(sps)
	cd = Dict(sps .=> colors[1:n])
	od = Dict()
	for (k, v) in d
		println(k)
		bend = pop!(bends)
		color = cd[specie(v)]
		style = "ultra thick, $color, bend $bend"
		od[_sort_tuple(k)] = style
	end
	return od
end

# ╔═╡ a932caba-ee96-11ea-040b-1530cd87a33e
sort(species(ab))

# ╔═╡ 4a0203a2-ee92-11ea-0aef-c105294c0613
edge_styles(g, rd1, ab)

# ╔═╡ 399c0a4e-ee81-11ea-04f5-613d49664145
tt = TikzGraphs.plot(g, node_styles=node_styles(g),
						edge_labels=edge_labels(g,rd1),
						edge_styles=edge_styles(g,rd1,ab),
						layout=Layouts.SpringElectrical(),
						options="scale=3")

# ╔═╡ b5e5766c-ee83-11ea-2708-09bdbc8cf736
begin
	TikzPictures.save(PDF("/Users/lyx/Downloads/testgraph"), tt)
	TikzPictures.save(SVG("/Users/lyx/Downloads/testgraph"), tt)
	TikzPictures.save(TEX("/Users/lyx/Downloads/testgraph"), tt)
end

# ╔═╡ 707ef074-ee85-11ea-2384-1be073e55fcb
plot_graph(ab; bends=[false, false, false, true, false, true])

# ╔═╡ 17ce563a-f007-11ea-1d65-4fded25ffb30
@show ab

# ╔═╡ 6b5665c0-f04a-11ea-3aa8-dbd665591deb
bcg = BlockCopolymerGraph(ab)

# ╔═╡ 889ec000-f04a-11ea-1fd1-a3c4e6dad652
ne(bcg.graph)

# ╔═╡ 8f157872-f04a-11ea-0a59-718fd6569e60
nv(bcg.graph)

# ╔═╡ e2f2a340-f04b-11ea-2d00-11def06a2166
chaintype(ab)

# ╔═╡ 717bbac0-f008-11ea-019d-69b77de3377e
diblock_chain()

# ╔═╡ 68799d18-f04a-11ea-0bc6-b70fcf410217
function AB3A3()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb0 = BranchPoint(:EB0)
	eb1 = BranchPoint(:EB1)
	eb2 = BranchPoint(:EB2)
	eb3 = BranchPoint(:EB3)
	A = PolymerBlock(:A, sA, 0.4, FreeEnd(:A), eb0)
	B1 = PolymerBlock(:B1, sB, 0.15, eb1, eb0)
	B2 = PolymerBlock(:B2, sB, 0.15, eb2, eb0)
	B3 = PolymerBlock(:B3, sB, 0.15, eb3, eb0)
	A1 = PolymerBlock(:A1, sA, 0.05, eb1, FreeEnd(:A1))
	A2 = PolymerBlock(:A2, sA, 0.05, eb2, FreeEnd(:A2))
	A3 = PolymerBlock(:A3, sA, 0.05, eb3, FreeEnd(:A3))
	return BlockCopolymer(:AB3A3, [A, B1, B2, B3, A1, A2, A3]) 
end

# ╔═╡ 47ba9282-f04d-11ea-0c27-c338d24c14a5
ab3a3 = AB3A3()

# ╔═╡ 52b167e2-f04d-11ea-0433-035113ba89f5
chaintype(ab3a3)

# ╔═╡ 65147104-f04d-11ea-3b9e-d127f95a74db
plot_graph(ab3a3; bends=[true, true, false, false, true, true, false])

# ╔═╡ Cell order:
# ╠═b6251b9a-ee70-11ea-2d4a-e56c3bcee76d
# ╠═80b5df48-ee76-11ea-1ce1-47c07d90a110
# ╠═88653aa4-ee76-11ea-1d32-59a8a26d2715
# ╠═27d3f3c4-ee85-11ea-31bd-8956d26c9c9b
# ╠═c090e90e-ee72-11ea-1b34-7b5c1ff84d24
# ╠═de11d6ae-ee73-11ea-380f-451eea158e9c
# ╠═0163aa4c-ee74-11ea-2f1f-b50463c30d24
# ╠═8707d6d8-ee73-11ea-1ce8-71fa50acb0bf
# ╠═df60ce48-f01c-11ea-0ccd-83995099b387
# ╠═821487a0-ee7a-11ea-2f7f-1de469ce28df
# ╠═0bde3d3a-ee7f-11ea-18cc-f3ee75c77ce1
# ╠═06ee1f2c-ee83-11ea-3f3a-adec2515cc8f
# ╠═2784cf78-ee83-11ea-363d-63a4fc2d0396
# ╠═39653f5a-ee83-11ea-0560-4d7ab9eaf280
# ╠═eebf4d96-ee83-11ea-37ed-01cb453eec53
# ╠═67e66760-ee91-11ea-14cc-07e697ba072c
# ╠═9fab2686-ee91-11ea-106d-b3c9ca2e4231
# ╠═040ebc18-ee84-11ea-163c-1fd9f5536e09
# ╠═a8bc6e66-ee84-11ea-0fcf-71f636966091
# ╠═07aaf17e-ee85-11ea-0723-677e9f5ae9e5
# ╠═a932caba-ee96-11ea-040b-1530cd87a33e
# ╠═4a0203a2-ee92-11ea-0aef-c105294c0613
# ╠═399c0a4e-ee81-11ea-04f5-613d49664145
# ╠═b5e5766c-ee83-11ea-2708-09bdbc8cf736
# ╠═707ef074-ee85-11ea-2384-1be073e55fcb
# ╠═17ce563a-f007-11ea-1d65-4fded25ffb30
# ╠═6b5665c0-f04a-11ea-3aa8-dbd665591deb
# ╠═889ec000-f04a-11ea-1fd1-a3c4e6dad652
# ╠═8f157872-f04a-11ea-0a59-718fd6569e60
# ╠═e2f2a340-f04b-11ea-2d00-11def06a2166
# ╠═717bbac0-f008-11ea-019d-69b77de3377e
# ╠═68799d18-f04a-11ea-0bc6-b70fcf410217
# ╠═47ba9282-f04d-11ea-0c27-c338d24c14a5
# ╠═52b167e2-f04d-11ea-0433-035113ba89f5
# ╠═65147104-f04d-11ea-3b9e-d127f95a74db
