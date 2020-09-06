### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ 50940824-edb5-11ea-1ae6-9bee2df52fd3
using TikzPictures

# ╔═╡ 7efbdc96-edb5-11ea-2167-59857741434e
using TikzGraphs

# ╔═╡ 86a79bec-edb5-11ea-3ad1-61586fa95455
using LightGraphs

# ╔═╡ 610245ea-edb5-11ea-2c63-072bdf57ae79
tp = TikzPicture("\\draw (0,0) -- (10,10);\n\\draw (10,0) -- (0,10);\n\\node at (5,5) {tikz \$\\sqrt{\\pi}\$};", options="scale=0.25", preamble="")

# ╔═╡ 8da31dfe-edb5-11ea-31f1-a99fbd18be84
g = DiGraph(4)

# ╔═╡ 97be7bd0-edb5-11ea-2c0f-2d942caa783b
add_edge!(g, 1, 2)

# ╔═╡ 9d57781c-edb5-11ea-3a3d-39da277ccf70
add_edge!(g, 2, 3)

# ╔═╡ a509c236-edb5-11ea-1977-0936cdb523af
TikzGraphs.plot(g)

# ╔═╡ b0d23e04-edb5-11ea-3504-8366dcd4115b
begin
	add_edge!(g, 3, 4)
	add_edge!(g, 1, 4)
	TikzGraphs.plot(g)
end

# ╔═╡ 0ae74ccc-edb6-11ea-1d07-39749f25f191
begin
	t = TikzGraphs.plot(g)
	TikzPictures.save(PDF("/Users/lyx/Downloads/testgraph"), t)
	TikzPictures.save(SVG("/Users/lyx/Downloads/testgraph"), t)
	TikzPictures.save(TEX("/Users/lyx/Downloads/testgraph"), t)
end

# ╔═╡ 20014eb8-edb7-11ea-356f-5f81ceb0c799
g2 = Graph(5)

# ╔═╡ 37645bcc-edb7-11ea-0ea8-3b75a7a4a187
begin
	add_edge!(g2, 1, 2)
	add_edge!(g2, 1, 3)
	add_edge!(g2, 1, 4)
	add_edge!(g2, 1, 5)
end

# ╔═╡ 7e131288-edb8-11ea-268d-3561f434460e
begin
	add_vertices!(g2, 3)
	add_edge!(g2, 3, 6)
	add_edge!(g2, 4, 7)
	add_edge!(g2, 5, 8)
end

# ╔═╡ 97cbe86c-edb8-11ea-3c2e-31ba7b2a07ea


# ╔═╡ 52bd8394-edb7-11ea-16a7-f74cb7141d60
TikzGraphs.plot(g2, Layouts.SimpleNecklace())

# ╔═╡ 7814f17c-edb7-11ea-2ed7-edfbdf3dd2a7
TikzGraphs.plot(g2, Layouts.Spring())

# ╔═╡ 9d04271e-edb7-11ea-1836-5f1fa8aebf88
TikzGraphs.plot(g2)

# ╔═╡ ab37a862-edb7-11ea-2bc9-6dcffb8330e3
TikzGraphs.plot(g2, ["∘","","∘","∘","∘","","",""],
	edge_styles=Dict((1,2)=>"red",
					 (1,3)=>"blue",
					 (1,4)=>"blue",
					 (1,5)=>"blue",
					 (3,6)=>"red",
					 (4,7)=>"red",
					 (5,8)=>"red"),
	layout=Layouts.Spring())

# ╔═╡ 5aaadb72-edb9-11ea-3672-ef34095d2be3
tt = TikzGraphs.plot(g2,
	node_styles=Dict(1=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 2=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 3=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 4=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 5=>"empty nodes, minimum size=0pt, inner sep=0pt",
				     6=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 7=>"empty nodes, minimum size=0pt, inner sep=0pt",
					 8=>"empty nodes, minimum size=0pt, inner sep=0pt"),
	edge_labels=Dict((1,2)=>"0.4",
					 (1,3)=>"0.15",
					 (1,4)=>"0.15",
					 (1,5)=>"0.15",
					 (3,6)=>"0.05",
					 (4,7)=>"0.05",
					 (5,8)=>"0.05"),
	edge_styles=Dict((1,2)=>"very thick, red, bend left",
					 (1,3)=>"very thick, blue, bend left",
					 (1,4)=>"very thick, blue, bend left",
					 (1,5)=>"very thick, blue, bend left",
					 (3,6)=>"very thick, red, bend right",
					 (4,7)=>"very thick, red, bend right",
					 (5,8)=>"very thick, red, bend right"),
	layout=Layouts.SpringElectrical(), options="scale=3")

# ╔═╡ 3ecb6320-edbb-11ea-2a5d-ef217f17f3c4
begin
	TikzPictures.save(PDF("/Users/lyx/Downloads/testgraph"), tt)
	TikzPictures.save(SVG("/Users/lyx/Downloads/testgraph"), tt)
	TikzPictures.save(TEX("/Users/lyx/Downloads/testgraph"), tt)
end

# ╔═╡ Cell order:
# ╠═50940824-edb5-11ea-1ae6-9bee2df52fd3
# ╠═610245ea-edb5-11ea-2c63-072bdf57ae79
# ╠═7efbdc96-edb5-11ea-2167-59857741434e
# ╠═86a79bec-edb5-11ea-3ad1-61586fa95455
# ╠═8da31dfe-edb5-11ea-31f1-a99fbd18be84
# ╠═97be7bd0-edb5-11ea-2c0f-2d942caa783b
# ╠═9d57781c-edb5-11ea-3a3d-39da277ccf70
# ╠═a509c236-edb5-11ea-1977-0936cdb523af
# ╠═b0d23e04-edb5-11ea-3504-8366dcd4115b
# ╠═0ae74ccc-edb6-11ea-1d07-39749f25f191
# ╠═20014eb8-edb7-11ea-356f-5f81ceb0c799
# ╠═37645bcc-edb7-11ea-0ea8-3b75a7a4a187
# ╠═7e131288-edb8-11ea-268d-3561f434460e
# ╠═97cbe86c-edb8-11ea-3c2e-31ba7b2a07ea
# ╠═52bd8394-edb7-11ea-16a7-f74cb7141d60
# ╠═7814f17c-edb7-11ea-2ed7-edfbdf3dd2a7
# ╠═9d04271e-edb7-11ea-1836-5f1fa8aebf88
# ╠═ab37a862-edb7-11ea-2bc9-6dcffb8330e3
# ╠═5aaadb72-edb9-11ea-3672-ef34095d2be3
# ╠═3ecb6320-edbb-11ea-2a5d-ef217f17f3c4
