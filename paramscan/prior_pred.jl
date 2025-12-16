using Pkg
Pkg.activate("..")

using Distributed

# Start workers
const WANT = 60
addprocs(max(1, min(WANT, Sys.CPU_THREADS - 1)))

@everywhere using Distributions, Agents, Statistics, DataFrames, CSV,
                  StatsBase, Random, Graphs, ProgressMeter, SharingSignals

@everywhere begin
    # Sim for prior predictive using the SAME summary definitions used in fitting
    function run_ppc_simulation(row::NamedTuple)
        tmax = 1000
        model = initialize_sharing_signals_ywb(
            l = row.l, C = row.C, γ = row.γ, beta = row.β,
            dens = row.dens, ref = row.ref, strength = row.strength,
            seed = Int(row.seed), total_ticks = tmax
        )
        for _ in 1:tmax
            step!(model)
        end

        avg_deg = mean(degree(model.sharenet))
        endows = [a.harvest for a in allagents(model)|>collect]
        
        return (
            share_avg_degree = avg_deg,
            share_median_indegree = median(indegree(model.sharenet)),
            share_median_outdegree = median(outdegree(model.sharenet)),
            share_avg_clust = mean([a.share_clust for a in allagents(model)]),
            share_coreness_avg = mean(core_number(model.sharenet)),
            share_var_indegree = var(indegree(model.sharenet)),
            share_iqr_indegree = iqr(indegree(model.sharenet)),
            share_mad_indegree = mad(indegree(model.sharenet)),
            share_var_outdegree = var(outdegree(model.sharenet)),
            share_iqr_outdegree = iqr(outdegree(model.sharenet)),
            share_mad_outdegree = mad(outdegree(model.sharenet)),
            share_zero_indegree_count = count(==(0), indegree(model.sharenet)),
            pagerank_harvest_corr = avg_deg > 0 ? model.pagerank_harvest_corr = corspearman([a.pagerank for a in allagents(model)|>collect], endows) : 0.0,
            gini_indegree = gini(indegree(model.sharenet)),
            reciprocity_dyad = reciprocity_igraph(model.sharenet),
            reciprocity_edge = reciprocity_igraph(model.sharenet, mode=:edge),
            com_endow_corr = model.dens > 0 ? corspearman([a.comdeg for a in allagents(model)|>collect], endows) : 0.0,
            c = model.C,
            l = model.l,
            gamma = model.γ,
            beta = model.beta,
            dens = model.dens,
            ref = model.ref,
            strength = model.strength
        )
    end
end

# --- LHS prior sampler (as you wrote)
@everywhere function lhs_sample_uniform(ns, rng)
    step = 1.0/ns
    pts  = step .* (collect(0:ns-1) .+ rand(rng, ns))
    shuffle!(rng, pts)
end
@everywhere lhs_sample(n, lo, hi, rng) = lo .+ (hi-lo).*lhs_sample_uniform(n, rng)
@everywhere lhs_sample_exponential(n, λ, rng) = -log.(lhs_sample_uniform(n, rng))./λ

@everywhere function sample_from_prior_lhs(N, rate, rng)
    l = lhs_sample_exponential(N, rate[1], rng)
    C = lhs_sample_exponential(N, rate[2], rng)
    γ = lhs_sample(N, 0, 1, rng)
    β = lhs_sample_exponential(N, rate[3], rng)
    dens = lhs_sample(N, 0, 1, rng)
    ref = lhs_sample(N, 0, 1, rng)
    strength = lhs_sample(N, 0, 1, rng)
    seeds = rand(rng, 1:10^8, N)
    # Return (l, C, γ, β, dens, ref, strength, seed)
    [(l[i], C[i], γ[i], β[i], dens[i], ref[i], strength[i], seeds[i]) for i in 1:N]
end

# --- Build the prior sample DataFrame
n_samples = 10_000
prior_tuples = sample_from_prior_lhs(n_samples, [1.0, 0.1, 1.0], Xoshiro(80085))
sampled_df = DataFrame(prior_tuples, [:l, :C, :γ, :β, :dens, :ref, :strength, :seed])

# Use NamedTuples to ship light rows to workers
rows_nt = collect(Tables.namedtupleiterator(sampled_df))

# --- Parallel prior predictive sims
chunk_size = max(1, nworkers())
nchunks = cld(n_samples, chunk_size)

p = Progress(n_samples, desc = "Prior predictive check simulations", dt = 1)
all_results = Vector{NamedTuple}()

for chunk_idx in 1:nchunks
    @info "Processing chunk $chunk_idx / $nchunks"
    s = 1 + (chunk_idx - 1)*chunk_size
    e = min(chunk_idx*chunk_size, n_samples)
    chunk_rows = rows_nt[s:e]

    # pmap across this chunk
    chunk_results = pmap(run_ppc_simulation, chunk_rows)
    append!(all_results, chunk_results)

    update!(p, e; force=true)  # absolute progress works fine
end

sim_netstats = DataFrame(all_results)

# Write results
output_ppc_path = "../data/prior_pred.csv"
CSV.write(output_ppc_path, sim_netstats)
@info "PPC simulations complete! Results saved to $output_ppc_path."