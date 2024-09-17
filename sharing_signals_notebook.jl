### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d62a86ec-7ba0-492c-b01e-93e87d514d28
begin
	using Revise
	using Pkg
	Pkg.activate(".")
	include("./src/SharingSignals_ABM.jl")
	using Agents
end

# ╔═╡ 9fe1740e-f211-4214-bd31-6113ebc90530
using StatsBase

# ╔═╡ 1a08fc7b-7cd7-4fa6-89ba-b0eecf23b161
using Random

# ╔═╡ 451f5b67-4c62-4e69-9b04-1d8a7af6e409
using Graphs

# ╔═╡ 0fa5d256-5fff-49dc-afd8-e29c88f70b5f
using Plots, GraphRecipes

# ╔═╡ b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
model = initialize_sharing_signals(B=[10, 2], proph=0.05, N=100, dens=0.2, C=1.0, b=0.1, beta=0.5, seed=rand(1:1000))

# ╔═╡ d3029894-faef-41d4-bcfe-9cedda7bd116
plot(0:25, social_reward.(0:25, beta=model.beta))

# ╔═╡ 623c1874-c0ed-4807-b002-a1a9f328a816
graphplot(model.comnet, curves=false)

# ╔═╡ 19903e9c-4ebf-4c6f-bf02-c59dc7f23f72
outdegree(model.sharenet)

# ╔═╡ 0df02e4d-b499-4f69-91a2-eb22f96e1fe7
run!(model, 10000, adata=[:payoff])

# ╔═╡ a65dcb32-4d99-4146-a7ab-4e43d7d32f0a
graphplot(model.sharenet, curves=false, node_weights=[a.payoff for a in allagents(model)|>collect])

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
	[a.B + indegree(model.sharenet, a.id)*model.C + a.payoff for a in allagents(model)|>collect]
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

# ╔═╡ 5d068d71-9569-44a3-82a8-a5fd15b48aa4
mean(outdegree(model.sharenet))

# ╔═╡ 80647b36-c52e-43db-806b-5f81266821dc
histogram([a.payoff for a in allagents(model)|>collect])

# ╔═╡ 48ca717e-90c0-4a9b-bb4f-832c869c96ec
mean([a.payoff for a in allagents(model)|>collect])

# ╔═╡ 4dd4a2b9-d919-4be1-97af-f070e8912db1
scatter(
	[a.B for a in allagents(model)|>collect],
	[a.payoff for a in allagents(model)|>collect]
)

# ╔═╡ b4bd93d0-1883-452f-8ec1-12b38dd2b61b
f(k, beta) = (1 - 1 / (1 - exp(beta)*(1 - model.dens)))*((1-(1-model.dens)^k)/model.dens) + (1 / (1 - exp(beta)*(1 - model.dens)))*(exp(beta)*(1 - exp(-beta * k)) / (exp(beta) - 1))

# ╔═╡ 2a19c151-6660-4042-8a9d-3aa0356df3a3
g(k) = (1 - (1 - model.dens)^k)*(1 / model.dens)

# ╔═╡ 83449d58-b4c5-4aff-9eb0-e600276ee156
h(k) = ( -(((1 - model.dens)^k - 1)*(model.dens - 1)) + k*model.dens ) / model.dens^2

# ╔═╡ 46f4f501-64ae-40e6-be13-c9f35c9fb344
f(mean(outdegree(model.sharenet)), model.beta)*(model.dens * model.N) - mean(outdegree(model.sharenet))*model.C

# ╔═╡ c9754d72-b756-4099-a793-f5ea58ee0ffc
g(mean(outdegree(model.sharenet)))*(model.dens * model.N) - mean(outdegree(model.sharenet))*model.C

# ╔═╡ 35cb145f-5442-40ea-a1a1-d9a71bf645bc
h(mean(outdegree(model.sharenet)))*(model.dens * model.N) - mean(outdegree(model.sharenet))*model.C

# ╔═╡ f45944fc-3051-4855-b381-7cd14e8b904b
neighbors(model.comnet, 3)

# ╔═╡ Cell order:
# ╠═d62a86ec-7ba0-492c-b01e-93e87d514d28
# ╠═9fe1740e-f211-4214-bd31-6113ebc90530
# ╠═1a08fc7b-7cd7-4fa6-89ba-b0eecf23b161
# ╠═451f5b67-4c62-4e69-9b04-1d8a7af6e409
# ╠═0fa5d256-5fff-49dc-afd8-e29c88f70b5f
# ╠═b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
# ╠═d3029894-faef-41d4-bcfe-9cedda7bd116
# ╠═623c1874-c0ed-4807-b002-a1a9f328a816
# ╠═19903e9c-4ebf-4c6f-bf02-c59dc7f23f72
# ╠═0df02e4d-b499-4f69-91a2-eb22f96e1fe7
# ╠═a65dcb32-4d99-4146-a7ab-4e43d7d32f0a
# ╠═c5ec000b-fc0b-4823-b672-430df33b0933
# ╠═8f4cda81-1255-4e27-b479-80f0f309f2a6
# ╠═1e6d1b71-e032-4acd-86a9-448707c74d0d
# ╠═04e5a2e0-5660-4075-9015-860fdfb9c2c0
# ╠═bf8b39bf-5f8a-470a-913f-84615bbcb185
# ╠═57e77f85-a8b5-4ef1-ade0-43bd53574a70
# ╠═5d068d71-9569-44a3-82a8-a5fd15b48aa4
# ╠═80647b36-c52e-43db-806b-5f81266821dc
# ╠═48ca717e-90c0-4a9b-bb4f-832c869c96ec
# ╠═4dd4a2b9-d919-4be1-97af-f070e8912db1
# ╠═b4bd93d0-1883-452f-8ec1-12b38dd2b61b
# ╠═2a19c151-6660-4042-8a9d-3aa0356df3a3
# ╠═83449d58-b4c5-4aff-9eb0-e600276ee156
# ╠═46f4f501-64ae-40e6-be13-c9f35c9fb344
# ╠═c9754d72-b756-4099-a793-f5ea58ee0ffc
# ╠═35cb145f-5442-40ea-a1a1-d9a71bf645bc
# ╠═f45944fc-3051-4855-b381-7cd14e8b904b
