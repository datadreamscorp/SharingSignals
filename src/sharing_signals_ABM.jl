using Statistics, StatsBase, Random, Distributions, Agents, Graphs, LinearAlgebra, SparseArrays

@agent struct Peep(NoSpaceAgent)
    harvest::Float64
    endow::Float64
    payoff::Float64
    total_payoff::Float64
    reps::Vector{Float64}
    signals::Vector{Bool}
    comneighbors::Vector{Int}
    comdeg::Int
    share_outneighbors::Vector{Int}
    outdeg::Int
    share_inneighbors::Vector{Int}
    indeg::Int
    com_clust::Float64
    share_clust::Float64
    pagerank::Float64
    core::Float64
end

Base.@kwdef mutable struct Parameters
    N::Int
    B::Float64
    sigma::Float64
    b::Float64
    C::Float64
    l::Float64
    h::Float64
    γ::Float64
    beta::Float64
    dens::Float64
    strength::Float64
    ref::Float64
    payoff_mode::Symbol = :exponential
    comnet::SimpleGraph
    rank::Bool
    sharenet::SimpleDiGraph
    tick::Int
    total_ticks::Int
    # DATA
    total_warns::Int = 0
    share_density::Float64 = 0.0
    com_density::Float64 = 0.0
    share_avg_degree::Float64 = 0.0
    share_median_indegree::Float64 = 0.0
    share_maximum_indegree::Float64 = 0.0
    share_median_outdegree::Float64 = 0.0
    com_avg_degree::Float64 = 0.0
    share_var_degree::Float64 = 0.0
    share_var_indegree::Float64 = 0.0
    share_iqr_indegree::Float64 = 0.0
    share_mad_indegree::Float64 = 0.0
    share_var_outdegree::Float64 = 0.0
    share_iqr_outdegree::Float64 = 0.0
    share_mad_outdegree::Float64 = 0.0
    share_avg_clust::Float64 = 0.0
    com_avg_clust::Float64 = 0.0
    share_giant_component::Int = 0
    com_giant_component::Int = 0
    share_diameter::Float64 = 0.0
    com_diameter::Float64 = 0.0
    share_avg_pathlength::Float64 = 0.0
    com_avg_pathlength::Float64 = 0.0
    com_degree_assort::Float64 = 0.0
    share_degree_assort::Float64 = 0.0
    com_wealth_assort::Float64 = 0.0
    share_wealth_assort::Float64 = 0.0
    #reciprocity
    share_reciprocity_dyad::Float64 = 0.0
    share_reciprocity_edge::Float64 = 0.0
    # correlation + zero-degree fields
    com_endow_deg_corr::Float64 = 0.0
    share_endow_indeg_corr::Float64 = 0.0
    share_endow_outdeg_corr::Float64 = 0.0
    com_zero_deg_count::Int = 0
    share_zero_outdegree_count::Int = 0
    share_zero_indegree_count::Int = 0
    com_endow_deg_corr_nonzero::Float64 = 0.0
    share_endow_deg_corr_nonzero::Float64 = 0.0
    # coreness
    com_coreness_avg::Float64 = 0.0
    com_coreness_max::Float64 = 0.0
    com_coreness_endow_corr::Float64 = 0.0
    share_coreness_avg::Float64 = 0.0
    share_coreness_max::Float64 = 0.0
    share_coreness_endow_corr::Float64 = 0.0
    pagerank_harvest_corr::Float64 = 0.0
    gini_indegree::Float64 = 0.0

end

@inline social_reward(x; beta = 0.1) = x > 0 ? (beta == 0 ? x : (1 - exp(-beta*x))/beta) : 0

function gini(v::AbstractVector{<:Real})
    n = length(v)
    n == 0 && return 0.0
    s = sort(float.(v))                # ascending
    total = sum(s)
    total == 0 && return 0.0
    coef = 2 / (n * total)
    g = coef * sum( j * s[j] for j in 1:n ) - (n + 1) / n
    return g                           # 0  ≤ g ≤ 1
end

function measure_lcc_diameter(g::SimpleGraph)
    if nv(g) <= 1 || ne(g) == 0
        return (nv(g), 0.0, 0.0)
    end
    comps = connected_components(g)
    # "giant = maximum(comps, by=length)" replaced by reduce():
    giant = reduce((x, y) -> length(x) >= length(y) ? x : y, comps)
    subg, mapping = induced_subgraph(g, giant)
    return (length(giant), diameter(subg), mean_distance(subg))
end

"""
    my_bfs_distances(g, start) -> Vector{Float64}

Returns a distance vector where `d[i]` is the number of edges in the
shortest path from `start` to node `i`. If a node is not reachable,
its distance will be Inf.
"""
function bfs_distances(g::AbstractGraph, start::Int)
    n = nv(g)
    dist = fill(Inf, n)
    dist[start] = 0
    visited = fill(false, n)
    visited[start] = true
    queue = [start]

    while !isempty(queue)
        v = popfirst!(queue)
        dv = dist[v]
        for w in neighbors(g, v)
            if !visited[w]
                visited[w] = true
                dist[w] = dv + 1
                push!(queue, w)
            end
        end
    end
    return dist
end

"""
    mmean_distance(g)

Compute the mean shortest-path distance between all connected pairs
of vertices in `g`, ignoring infinite distances (i.e. pairs that are
not in the same connected component).
"""
function mean_distance(g::AbstractGraph)
    n = nv(g)
    totaldist = 0.0
    count = 0
    for v in 1:n
        distv = bfs_distances(g, v)
        for w in 1:n
            d = distv[w]
            if isfinite(d) && w != v
                totaldist += d
                count += 1
            end
        end
    end
    return count == 0 ? 0.0 : totaldist / count
end

function measure_lcc_diameter(g::SimpleDiGraph)
    ug = SimpleGraph(g)
    return measure_lcc_diameter(ug)
end

function assortativity_by_attribute_undirected(g::SimpleGraph, x::AbstractVector{<:Real})
    E = edges(g)
    if isempty(E)
        return 0.0
    end
    xv = [x[src(e)] for e in E]
    yv = [x[dst(e)] for e in E]
    return cor(xv, yv)
end


@inline function assortativity_by_attribute_directed(g::SimpleDiGraph, x::AbstractVector{<:Real})
    ug = SimpleGraph(g)
    return assortativity_by_attribute_undirected(ug, x)
end

@inline function degree_assortativity_com(model)
    degs = [degree(model.comnet, i) for i in 1:model.N]
    return assortativity_by_attribute_undirected(model.comnet, degs)
end

@inline function degree_assortativity_share(model)
    ug = SimpleGraph(model.sharenet)
    degs = [degree(ug, i) for i in 1:model.N]
    return assortativity_by_attribute_undirected(ug, degs)
end

@inline function wealth_assortativity_com(model)
    wealths = [model[i].harvest for i in 1:model.N]
    return assortativity_by_attribute_undirected(model.comnet, wealths)
end

@inline function wealth_assortativity_share(model)
    wealths = [model[i].harvest for i in 1:model.N]
    return assortativity_by_attribute_directed(model.sharenet, wealths)
end

function reciprocity_igraph(g; mode::Symbol = :dyad)
    # sparse adjacency matrix of g (out-edges => rows are sources, cols are targets)
    A  = adjacency_matrix(g; dir = :out)
    AT = transpose(A)

    if mode === :edge
        # edge-based: sum(A .* Aᵀ) counts EACH orientation of a mutual dyad
        recip_edges = sum(A .* AT)
        return nnz(A) == 0 ? NaN : recip_edges / nnz(A)

    elseif mode === :dyad
        # dyad-based: each mutual pair counted ONCE
        mutual_dyads = sum(A .* AT) ÷ 2 # divide by 2 because A.*Aᵀ hits both (i,j) and (j,i)
        # any-edge dyads: look at strict lower triangle (i>j) of A+AT
        present_dyads = count(!iszero, tril(A + AT, -1))
        return present_dyads == 0 ? NaN : mutual_dyads / present_dyads

    else
        throw(ArgumentError("mode must be :dyad or :edge"))
    end
end

function gather_network_stats!(model)
    n = model.N
    
    # 1) DENSITY + AVERAGE DEGREE
    edges_com = ne(model.comnet)
    edges_share = ne(model.sharenet)
    
    com_net_density = (n > 1) ? (2 * edges_com    / (n*(n-1))) : 0.0
    share_net_density = (n > 1) ? (   edges_share   / (n*(n-1))) : 0.0
    
    model.com_density = com_net_density
    model.share_density = share_net_density
    model.com_avg_degree = mean(degree(model.comnet))
    model.share_avg_degree = mean(degree(model.sharenet))
    model.share_median_indegree = median(indegree(model.sharenet))
    model.share_maximum_indegree = maximum(indegree(model.sharenet))
    model.share_median_outdegree = median(outdegree(model.sharenet))
    model.share_var_degree = var(degree(model.sharenet))
    model.share_var_indegree = var(indegree(model.sharenet))
    model.share_iqr_indegree = iqr(indegree(model.sharenet))
    model.share_mad_indegree = mad(indegree(model.sharenet))
    model.share_var_outdegree = var(outdegree(model.sharenet))
    model.share_iqr_outdegree = iqr(outdegree(model.sharenet))
    model.share_mad_outdegree = mad(outdegree(model.sharenet))

    # AVERAGE LOCAL CLUSTERING
    agent_clust_com = [a.com_clust   for a in allagents(model)]
    agent_clust_share = [a.share_clust for a in allagents(model)]
    model.com_avg_clust = mean(agent_clust_com)
    model.share_avg_clust = mean(agent_clust_share)

    # ASSORTATIVITY
    model.com_degree_assort = degree_assortativity_com(model)
    model.share_degree_assort = degree_assortativity_share(model)
    model.com_wealth_assort = wealth_assortativity_com(model)
    model.share_wealth_assort = wealth_assortativity_share(model)

    # reciprocity
    model.share_reciprocity_dyad = reciprocity_igraph(model.sharenet)
    model.share_reciprocity_edge = reciprocity_igraph(model.sharenet, mode=:edge)

    # ENDOWMENT–DEGREE CORRELATIONS & ZERO-DEGREE COUNTS
    # collect degrees for each network
    com_degs   = [degree(model.comnet, i) for i in 1:n]
    share_indegs = [model[i].indeg for i in 1:n]
    share_outdegs = [model[i].outdeg for i in 1:n]
    endows = [model[i].harvest for i in 1:n]

    # correlation of endowment with communication-degree
    model.com_endow_deg_corr = corspearman(com_degs, endows)

    all_same(v) = isempty(v) || all(==(v[1]), v)

    # correlation of endowment with total sharing-degree
    model.share_endow_indeg_corr = all_same(share_indegs) ? 0.0 : corspearman(share_indegs, endows)
    model.share_endow_outdeg_corr = all_same(share_outdegs) ? 0.0 : corspearman(share_outdegs, endows)

    # zero-degree counts
    model.com_zero_deg_count = count(==(0), com_degs)
    model.share_zero_outdegree_count = count(==(0), [model[i].outdeg for i in 1:n])
    model.share_zero_indegree_count = count(==(0), [model[i].indeg  for i in 1:n])

    # K-CORE STRUCTURE
    # communication network coreness
    com_corevals = core_number(model.comnet)
    model.com_coreness_avg = mean(com_corevals)
    model.com_coreness_max = maximum(com_corevals)
    
    # correlation of coreness with endowment
    endows = [a.harvest for a in allagents(model)|>collect]
    model.com_coreness_endow_corr = corspearman(com_corevals, endows)

    # sharing network coreness
    ug_share = model.sharenet  # convert to undirected
    share_corevals = core_number(ug_share)
    model.share_coreness_avg = mean(share_corevals)
    model.share_coreness_max = maximum(share_corevals)

    # correlation of coreness and pagerank with harvest
    if model.share_avg_degree > 0
        model.share_coreness_endow_corr = corspearman(share_corevals, endows)
        model.pagerank_harvest_corr = corspearman([a.pagerank for a in allagents(model)|>collect], endows)
        model.gini_indegree = gini(indegree(model.sharenet))
    end
    
end

function payoff_exponential!(a, model)
    a.payoff = 0.0
    total_neighbors = Set{Int}()
    for n in a.share_outneighbors
        for t in model[n].comneighbors
            if t ∉ total_neighbors
                a.payoff += model[t].reps[a.id]
                push!(total_neighbors, t)
            end
        end
    end
    cost = (model.C*length(a.share_outneighbors)^model.l) / ( a.harvest^model.γ )
    # subtract the cost term
    a.payoff = max(0.0, a.payoff - cost)
end


function payoff_log!(a, model)
    # social reward
    soc_reward = 0.0
    total_neighbors = Set{Int}()
    for n in a.share_outneighbors
        for t in model[n].comneighbors
            if t ∉ total_neighbors
                soc_reward += model[t].reps[a.id]
                push!(total_neighbors, t)
            end
        end
    end

    cost = (length(a.share_outneighbors)^model.l) / (a.harvest)^model.γ
    a.payoff = max(0.0, soc_reward - cost)
end


function payoff!(a, model)
    if model.payoff_mode == :exponential
        payoff_exponential!(a, model)
    elseif model.payoff_mode == :log
        payoff_log!(a, model)
    else
        error("Unknown payoff_mode: $(model.payoff_mode). Use :exponential or :log.")
    end
end


function potential_payoff_exponential(a, j, model)
    pay = 0.0
    total_neighbors = Set{Int}()
    for n in a.share_outneighbors
        for t in model[n].comneighbors
            if t ∉ j.comneighbors && t ∉ total_neighbors
                pay += model[t].reps[a.id]
                push!(total_neighbors, t)
            end
        end
    end

    newrep = 0.0
    for t in j.comneighbors
        ssum = 0
        for m in model[t].comneighbors
            ssum += model[m].signals[a.id] ? 1 : 0
        end
        newrep = social_reward(ssum + 1, beta=model.beta)
        pay += newrep
    end

    cost = (model.C*(length(a.share_outneighbors) + 1)^model.l) / ( a.harvest^model.γ )

    return max(0.0, pay - cost)
end

function potential_payoff_log(a, j, model)
    pay = 0.0
    total_neighbors = Set{Int}()
    for n in a.share_outneighbors
        for t in model[n].comneighbors
            if t ∉ j.comneighbors && t ∉ total_neighbors
                pay += model[t].reps[a.id]
                push!(total_neighbors, t)
            end
        end
    end

    # newrep from j.comneighbors
    for t in j.comneighbors
        ssum = 0
        for m in model[t].comneighbors
            ssum += model[m].signals[a.id] ? 1 : 0
        end
        local_rep = social_reward(ssum + 1, beta=model.beta)
        pay += local_rep
    end

    # cost same as log payoff approach
    cost = ((length(a.share_outneighbors)+1)^model.l) / (a.harvest)^model.γ

    return max(0.0, pay - cost)
end

function potential_payoff(a, j, model)
    if model.payoff_mode == :exponential
        return potential_payoff_exponential(a, j, model)
    elseif model.payoff_mode == :log
        return potential_payoff_log(a, j, model)
    else
        error("Unknown payoff_mode $(model.payoff_mode). Use :exponential or :log.")
    end
end

function change_impression!(a, chosen, model)
    for n in chosen.comneighbors
        # sum of signals from model[n].comneighbors
        ssum = 0
        for m in model[n].comneighbors
            ssum += model[m].signals[a.id] ? 1 : 0
        end
        model[n].reps[a.id] = social_reward(ssum, beta=model.beta)
    end
end

function connect!(model)
    # Shuffle the agents in place
    ags = collect(allagents(model))
    shuffle!(abmrng(model), ags)

    all_ids = collect(1:model.N)  # or allids(model), if that’s the same

    for a in ags
        k = a.outdeg

        # check if we can afford to add another out-link
        if a.harvest > 0
            # build candidate list
            taken = Set([a.id])
            union!(taken, a.share_outneighbors)
            candidates = Int[]
            for c in all_ids
                if c ∉ taken
                    push!(candidates, c)
                end
            end

            if !isempty(candidates)
                chosen_agent = if model.rank
                    # sort by comdeg and pick the last
                    # sort(...) can be expensive, so just do a simple maximum:
                    max_agent = nothing
                    max_deg = -1
                    for cid in candidates
                        degval = model[cid].comdeg
                        if degval > max_deg
                            max_deg = degval
                            max_agent = cid
                        end
                    end
                    model[max_agent]
                else
                    # pick random from candidates
                    model[rand(abmrng(model), candidates)]
                end

                # if potential payoff is better than current total_payoff, link them
                if a.payoff < potential_payoff(a, chosen_agent, model)
                    add_edge!(model.sharenet, a.id, chosen_agent.id)
                    push!(a.share_outneighbors, chosen_agent.id)
                    a.outdeg += 1
                    push!(chosen_agent.share_inneighbors, a.id)
                    chosen_agent.indeg += 1
                    chosen_agent.signals[a.id] = true
                    change_impression!(a, chosen_agent, model)
                    payoff!(a, model)
                    payoff!(chosen_agent, model)
                end
            end
        end
    end

    pr   = pagerank(model.sharenet)
    core = core_number(model.sharenet)

    for a in ags
        a.com_clust   = local_clustering_coefficient(model.comnet, a.id)
        a.share_clust = local_clustering_coefficient(model.sharenet, a.id)
        a.pagerank = pr[a.id]
        a.core = core[a.id]
    end
    
    model.tick += 1
    #model.tick == model.total_ticks && gather_network_stats!(model)
end

function initialize_sharing_signals_ywb(;
    N = 100,
    B = 1.0,
    sigma = 3.0,
    b = 1.0,
    C = 1.0,
    l = 1.0,
    h = 1.0,
    γ = 0.0,
    beta = 0.01,
    dens = 0.05,
    net_type = "random",
    payoff_mode = :exponential,
    strength = 0.0,
    ref = 0.0,
    rank = false,
    seed = 75648,
    total_ticks = 1000,
    steps = 1000,
)
    rng = Xoshiro(seed)

    endist = [
        782950, 445110, 0, 0, 9675670, 741850, 288870, 1043250, 0, 0, 0, 169630, 
        1081550, 858330, 0, 726570, 1276950, 402300, 832640, 318450, 697950, 0, 
        501310, 1865550, 56700, 176470, 698000, 0, 242690, 1789980, 3485440, 0, 0, 0, 
        364680, 0, 28100, 28100, 523950, 148370, 202260, 895680, 7485180, 0, 2425760, 
        32180, 876380, 853100, 1378660, 4458230, 0, 1514280, 0, 2299130, 0, 1011600, 
        2761740, 0, 0, 2077530, 0, 3995820, 521480, 2996530, 0, 522620, 0, 0, 1814230, 
        328920, 76370, 0, 0, 151180, 0, 1439370, 1169210, 0, 0, 148370, 148370, 
        180870, 148370, 0, 0, 296740, 703120, 472440, 153990, 156590, 1294740, 614740, 
        148370, 0, 1243450, 0, 1289760, 4108220, 535100, 0, 0, 0, 1278650, 895740, 
        296740, 0, 0, 0, 599100, 84300
    ] ./ 2500

    endowments = isempty(endist) ? rand(rng, LogNormal(log(B), log(sigma)), N).^(h) : endist.^h
    n = isempty(endist) ? N : length(endist)

    # create net
    if net_type == "random"
        if dens < 1
            net, correlation, total_warns = run_network_simulation(n, dens, log.(1 .+ endist), steps, strength, ref, rng=rng)
        else
            net = complete_graph(n)
            correlation = 0.0
            total_warns = 0
        end
    end

    properties = Parameters(
        N = n, B = B, sigma = sigma, b = b, C = C, l = l, h = h, γ = γ, beta = beta,
        dens = dens, strength = strength, ref = ref, comnet = net, rank = rank,
        sharenet = SimpleDiGraph(n),
        tick = 0, total_ticks = total_ticks, total_warns = total_warns
    )

    model = StandardABM(
        Peep, 
        nothing;
        properties = properties,
        model_step! = connect!,
        rng = rng
    )

    model.N = length(endist)

    for a in 1:n
        agent = Peep(
            a, # id
            endist[a], # harvest (just reusing input)
            endowments[a], # endow
            0.0, # payoff
            endowments[a], # total_payoff
            fill(0.0, n), # reps
            fill(false, n), # signals
            neighbors(net, a), # comneighbors
            degree(net, a), # comdeg
            Int[], # share_outneighbors
            0, # outdeg
            Int[], # share_inneighbors
            0, # indeg
            0.0, # com_clust
            0.0, # share_clust,
            0.0, # pagerank
            0.0 # core
        )
        new_a = add_agent!(agent, model)
        payoff!(new_a, model)
    end

    model.tick += 1

    return model
end

###
# network generation routines
###
@inline function transform_endowments(endowments::Vector{Float64}, reference::Float64)
    norm_ref = reference * maximum(endowments)
    return exp.(-abs.(endowments .- norm_ref))
end

function choose_pair_add(w, g, S, rng)
    rng = isnothing(rng) ? default_rng() : rng

    ne = [(i,j) for i in 1:nv(g) for j in i+1:nv(g) if !has_edge(g,i,j)]
    isempty(ne) && return nothing

    # vectorised weights
    wi = w[first.(ne)]
    wj = w[last.(ne)]
    probs = (wi .* wj).^S

    return sample(rng, ne, Weights(probs))
end

function choose_edge_remove(w, g, S, rng)
    rng = isnothing(rng) ? default_rng() : rng

    es = collect(edges(g))
    isempty(es) && return nothing

    probs = ((1 .- w[src.(es)]) .* (1 .- w[dst.(es)])).^S
    return sample(rng, es, Weights(probs))
end

function evolve_network(g::Graph, endowments::Vector{Float64}, steps::Int, S::Float64, reference::Float64; rng=nothing)
    transformed_endowments = transform_endowments(endowments, reference)
    total_warns = 0
    for _ in 1:steps
        node_pair = choose_pair_add(transformed_endowments, g, S, rng)
        if node_pair !== nothing
            node1, node2 = node_pair
            add_edge!(g, node1, node2)
        else
            # @warn "choose_pair_add exceeded max_attempts; skipping edge addition."
            total_warns += 1
        end
        edge_to_remove = choose_edge_remove(transformed_endowments, g, S, rng)
        if node_pair !== nothing && edge_to_remove !== nothing
            rem_edge!(g, edge_to_remove)
        end
    end
    return g, total_warns
end

function run_network_simulation(N::Int, u::Float64, endowments::Vector{Float64}, steps::Int,
                                S::Float64, reference::Float64; rng=nothing)
    g = erdos_renyi(N, u, rng=rng)
    g, total_warns = evolve_network(g, endowments, steps, S, reference, rng=rng)
    degrees = degree(g)
    correlation = cor(degrees, endowments)
    return g, correlation, total_warns
end

#=
function initialize_sharing_signals(;
    N = 100,
    B = 1.0,
    sigma = 3.0,
    b = 1.0,
    C = 1.0,
    l = 1.0,
    h = 1.0,
    γ = 0.0,
    beta = 0.01,
    dens = 0.05,
    net_type = "random",
    payoff_mode = :exponential,
    strength = 0.0,
    ref = 0.0,
    rank = false,
    endist = Vector{Float64}(),
    seed = 75648,
    total_ticks = 1000,
    steps = 1000,
)

    rng = Xoshiro(seed)
    endowments = isempty(endist) ? rand(rng, LogNormal(log(B), log(sigma)), N).^(h) : endist.^h
    n = isempty(endist) ? N : length(endist)

    # create net
    if net_type == "random"
        if dens < 1
            net, correlation = run_network_simulation(n, dens, log(1 .+ endist), steps, strength, ref, rng=rng)
        else
            net = complete_graph(n)
        end
    end

    properties = Parameters(
        N = n, B = B, sigma = sigma, b = b, C = C, l = l, h = h, γ = γ, beta = beta,
        dens = dens, strength = strength, ref = ref, comnet = net, rank = rank,
        sharenet = SimpleDiGraph(n),
        tick = 0, total_ticks = total_ticks,
    )

    model = StandardABM(
        Peep, 
        nothing;
        properties = properties,
        model_step! = connect!,
        rng = rng
    )

    for a in 1:n
        agent = Peep(
            a, # id
            endist[a], # harvest (just reusing input)
            endowments[a], # endow
            0.0, # payoff
            endowments[a], # total_payoff
            fill(0.0, n), # reps
            fill(false, n), # signals
            neighbors(net, a), # comneighbors
            degree(net, a), # comdeg
            Int[], # share_outneighbors
            0, # outdeg
            Int[], # share_inneighbors
            0, # indeg
            0.0, # com_clust
            0.0, # share_clust,
            0.0, # pagerank
            0.0 # core
        )
        new_a = add_agent!(agent, model)
        payoff!(new_a, model)
    end

    return model
end
=#