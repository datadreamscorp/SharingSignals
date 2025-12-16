using Pkg
Pkg.activate("..")

using Distributed

# workers
WANT = 80
addprocs(max(1, min(WANT, Sys.CPU_THREADS - 1)))

@everywhere using Distributions, Agents, Statistics, DataFrames, CSV,
                  StatsBase, Random, Graphs, ProgressMeter, SharingSignals

@everywhere begin
    # Simulation for PPC using the SAME summary definitions you used to fit
    function run_ppc_simulation(row::NamedTuple, seed::Int)
        tmax = 1000
        model = initialize_sharing_signals_ywb(
            l = row.l, C = row.C, γ = row.γ, beta = row.beta,
            dens = row.dens, ref = row.ref, strength = row.strength,
            seed = seed, total_ticks = tmax
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
            pagerank_harvest_corr = avg_deg > 0 ? corspearman([a.pagerank for a in allagents(model)|>collect], endows) : 0.0,
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

# --- draw posterior samples properly (final round + weights)
abc_path = "../data/var_fit_pooled.csv"
abc_df   = CSV.read(abc_path, DataFrame)

last_r   = maximum(abc_df.round)
post_end = abc_df[abc_df.round .== last_r, :]

w = collect(post_end.weight)
w ./= sum(w)

rng = Xoshiro(80085)
n_samples = 10_000
sample_indices = sample(rng, 1:nrow(post_end), Weights(w), n_samples; replace=true)
sampled = post_end[sample_indices, [:l, :C, :γ, :beta, :dens, :ref, :strength]]

# use NamedTuples to ship light rows to workers
rows_nt = collect(Tables.namedtupleiterator(sampled))

# --- parallel PPC
chunk_size = max(1, nworkers())
nchunks = cld(n_samples, chunk_size)

p = Progress(n_samples; desc="Posterior predictive check simulations", dt=1)
all_results = Vector{NamedTuple}()

for chunk_idx in 1:nchunks
    s = 1 + (chunk_idx - 1) * chunk_size
    e = min(chunk_idx * chunk_size, n_samples)
    chunk_rows = rows_nt[s:e]
    seeds = rand(rng, 1:10^8, length(chunk_rows))  # match length precisely

    # pmap zips the two iterables
    chunk_results = pmap(run_ppc_simulation, chunk_rows, seeds)
    append!(all_results, chunk_results)

    update!(p, e; force=true)
end

sim_netstats = DataFrame(all_results)
output_ppc_path = "../data/post_pred_var.csv"
CSV.write(output_ppc_path, sim_netstats)
@info "PPC simulations complete! Results saved to $output_ppc_path."