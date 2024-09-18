### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 887a17bf-76fd-445b-a82d-ff5c5e45d3f9
begin
    using Pkg
    Pkg.activate(".")
    using Revise
	using SharingSignals
    using Agents, Graphs, Statistics
	using Plots, GraphRecipes
end

# ╔═╡ b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
model = initialize_sharing_signals(
	N=50, 
	B=1.0, 
	sigma=10.0, 
	b=0.5,
	C=0.5,
	beta=1.0,
	dens=0.1,
	seed=rand(1:1000)
)

# ╔═╡ b4bd93d0-1883-452f-8ec1-12b38dd2b61b
f(k, beta) = (1 - 1 / (1 - exp(beta)*(1 - model.dens)))*((1-(1-model.dens)^k)/model.dens) + (1 / (1 - exp(beta)*(1 - model.dens)))*(exp(beta)*(1 - exp(-beta * k)) / (exp(beta) - 1))

# ╔═╡ 2a19c151-6660-4042-8a9d-3aa0356df3a3
g(k) = (1 - (1 - model.dens)^k)*(1 / model.dens)

# ╔═╡ 83449d58-b4c5-4aff-9eb0-e600276ee156
h(k) = ( -(((1 - model.dens)^k - 1)*(model.dens - 1)) + k*model.dens ) / model.dens^2

# ╔═╡ 0df02e4d-b499-4f69-91a2-eb22f96e1fe7
begin
	run!(model, 1000, adata=[:payoff])
	plot(
		graphplot(model.comnet, curves=false, title="communication"),
		graphplot(model.sharenet, curves=false, node_weights=[a.payoff for a in allagents(model)|>collect], title="sharing"),
		layout=(1,2)
	)
end

# ╔═╡ 37543d3a-de17-4372-9077-73ed7ac1f889
histogram([a.B for a in allagents(model)|>collect])

# ╔═╡ c5ec000b-fc0b-4823-b672-430df33b0933
scatter(
	degree(model.comnet),
	indegree(model.sharenet)
)

# ╔═╡ 8f4cda81-1255-4e27-b479-80f0f309f2a6
cor(
	degree(model.comnet),
	indegree(model.sharenet)
)

# ╔═╡ 1e6d1b71-e032-4acd-86a9-448707c74d0d
scatter(
	[outdegree(model.sharenet, a.id) for a in allagents(model)|>collect],
	[a.total_payoff for a in allagents(model)|>collect]
)

# ╔═╡ 04e5a2e0-5660-4075-9015-860fdfb9c2c0
scatter(
	[a.B for a in allagents(model)|>collect],
	[outdegree(model.sharenet, a.id) for a in allagents(model)|>collect]
)

# ╔═╡ bf8b39bf-5f8a-470a-913f-84615bbcb185
histogram(model.sharenet|>indegree)

# ╔═╡ 57e77f85-a8b5-4ef1-ade0-43bd53574a70
histogram(outdegree(model.sharenet))

# ╔═╡ 80647b36-c52e-43db-806b-5f81266821dc
histogram([a.payoff for a in allagents(model)|>collect])

# ╔═╡ 48ca717e-90c0-4a9b-bb4f-832c869c96ec
mean([a.payoff for a in allagents(model)|>collect])

# ╔═╡ 46f4f501-64ae-40e6-be13-c9f35c9fb344
f(mean(outdegree(model.sharenet)), model.beta)*(model.dens * model.N)*model.b - mean(outdegree(model.sharenet))*model.C

# ╔═╡ 4dd4a2b9-d919-4be1-97af-f070e8912db1
scatter(
	[a.B for a in allagents(model)|>collect],
	[a.payoff for a in allagents(model)|>collect]
)

# ╔═╡ Cell order:
# ╠═887a17bf-76fd-445b-a82d-ff5c5e45d3f9
# ╠═b4bd93d0-1883-452f-8ec1-12b38dd2b61b
# ╠═2a19c151-6660-4042-8a9d-3aa0356df3a3
# ╠═83449d58-b4c5-4aff-9eb0-e600276ee156
# ╠═b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
# ╟─0df02e4d-b499-4f69-91a2-eb22f96e1fe7
# ╠═37543d3a-de17-4372-9077-73ed7ac1f889
# ╠═c5ec000b-fc0b-4823-b672-430df33b0933
# ╠═8f4cda81-1255-4e27-b479-80f0f309f2a6
# ╠═1e6d1b71-e032-4acd-86a9-448707c74d0d
# ╠═04e5a2e0-5660-4075-9015-860fdfb9c2c0
# ╠═bf8b39bf-5f8a-470a-913f-84615bbcb185
# ╠═57e77f85-a8b5-4ef1-ade0-43bd53574a70
# ╠═80647b36-c52e-43db-806b-5f81266821dc
# ╠═48ca717e-90c0-4a9b-bb4f-832c869c96ec
# ╠═46f4f501-64ae-40e6-be13-c9f35c9fb344
# ╠═4dd4a2b9-d919-4be1-97af-f070e8912db1
