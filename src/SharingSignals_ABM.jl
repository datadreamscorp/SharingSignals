using StatsBase, Random, Distributions, Agents, Graphs

@agent struct Peep(NoSpaceAgent)
    ###
    B::Float64 #endowment
    payoff::Float64
    reps::Vector{Float64}
    signals::Vector{Bool}
end


Base.@kwdef mutable struct Parameters
    #model parameters
    N::Int64
    B::Vector{Float64}
    b::Float64
    C::Float64
    beta::Float64
    dens::Float64
    proph::Float64
    comnet::SimpleGraph
    sharenet::SimpleDiGraph
    #data
    K_dist::Vector{Float64}
	Kbar::Float64
	V_dist::Vector{Float64}
	Vbar::Float64
    tick::Int64
end


social_reward(x; beta = 0.1) = x > 0 ? sum( [exp(-(n-1)*beta) for n in 1:x] ) : 0

function payoff!(a, model)
    a.payoff = 0 #+ indegree(model.sharenet, a.id)*model.C
    total_neighbors = []
    for n in outneighbors(model.sharenet, a.id)
        current_neighbors = neighbors(model.comnet, n)
        for t in filter(x -> x ∉ total_neighbors, current_neighbors)
            a.payoff += model[t].reps[a.id]*model.b
        end
        total_neighbors = vcat(total_neighbors, current_neighbors)
        a.payoff -= model.C
    end
end

function potential_payoff(a, j, model)
    pay = 0 #+ indegree(model.sharenet, a.id)*model.C
    receivers = outneighbors(model.sharenet, a.id)
    total_neighbors = []
    for n in receivers
        current_neighbors = neighbors(model.comnet, n)
        for t in filter(x -> x ∉ vcat(neighbors(model.comnet, j.id), total_neighbors), neighbors(model.comnet, n))
            pay += model[t].reps[a.id]*model.b
        end
        total_neighbors = vcat(total_neighbors, current_neighbors)
    end
    pay -= length(receivers)*model.C

    newrep = 0
    for t in neighbors(model.comnet, j.id)
        newrep = social_reward(
            sum(
                [
                    model[m].signals[a.id]        
                    for m in neighbors(model.comnet, t)
                ]
            ) + 1,
            beta = model.beta
            )
        pay += newrep*model.b
    end    
    return pay - model.C
end

function change_impression!(a, j, model)
    for n in neighbors(model.comnet, j.id)
        model[n].reps[a.id] = social_reward(
            sum(
                [
                    model[m].signals[a.id]
                    for m in neighbors(model.comnet, n)
                ]
            ),
            beta = model.beta
        )
    end
end

function connect!(model)
    for a in shuffle( abmrng(model), allagents(model)|>collect )

        k = outdegree(model.sharenet)[a.id]
        if a.B - (k+1)*model.C > 0
            taken = vcat([a.id], outneighbors(model.sharenet, a.id))
            candidates = filter(p -> p ∉ taken, allids(model)|>collect)
            if length(candidates) > 0
                
                chosen = model[rand(abmrng(model), candidates)]

                if a.payoff < potential_payoff(a, chosen, model)
                    add_edge!(model.sharenet, a.id, chosen.id)
                    chosen.signals[a.id] = true
                    change_impression!(a, chosen, model)
                    payoff!(a, model)
                    payoff!(chosen, model)
                end
            end
        end
        
    end
end

function initialize_sharing_signals(;
	N = 50,
    B = [5, 1],
    b = 1,
    C = 0.1,
    beta = 2,
    dens = 0.5,
    proph = 0.1,
    seed = 75648,
)
	rng = Xoshiro(seed)
	
	properties = Parameters(
		#model parameters
	    N,
        B,
        b,
        C,
        beta,
        dens,
        proph,
        dens < 1 ? erdos_renyi(N, dens, seed=seed) : complete_graph(N),
        SimpleDiGraph(N),
        Vector{Float64}(),
        0.0,
        Vector{Float64}(),
        0.0,
        0
	)

	model = StandardABM( 
		Peep, 
		nothing;
		properties = properties,
        model_step! = connect!,
		rng = rng
	)

	for a in 1:N
		
		group = rand(abmrng(model)) < model.proph ? 1 : 2

		agent = Peep( 
			a,
            B[group],
			0.0,
            repeat([0.0], N),
            repeat([false], N)
		)
        new_a = add_agent!(agent, model)
		
        payoff!(new_a, model)
	end
	
	return model
		
end