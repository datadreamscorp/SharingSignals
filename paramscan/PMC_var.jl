using Pkg
Pkg.activate("..")

using Distributed

WANT = 80  # number of worker processes
#OUTFILE = "../data/var_fit1.csv"
#OUTFILE = "../data/var_fit2.csv"
#OUTFILE = "../data/var_fit3.csv"
OUTFILE = "../data/var_fit4.csv"
#OUTFILE = "../data/var_fit5.csv"

addprocs(max(1, min(WANT, Sys.CPU_THREADS - 1)))

# make packages and constants available on all workers
@everywhere using SharingSignals, Distributions, Agents, Random, Graphs, 
                  StatsBase, LinearAlgebra, ProgressMeter, CSV, DataFrames

@everywhere const R_REPLICATES = 5

# define our smoothing kernel bandwidth struct
@everywhere mutable struct KernelBandwidth
    σ_l::Float64
    σ_C::Float64
    σ_γ::Float64
    σ_beta::Float64
    σ_dens::Float64
    σ_ref::Float64
    σ_strength::Float64
end

# one‐run simulation and summary
@everywhere function run_model_expt(params::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,Int64})
    l, C, γ, β, dens, ref, strength, seed = params
    model = initialize_sharing_signals_ywb(
        l = l, C = C, γ = γ, beta = β,
        dens = dens, ref = ref, strength = strength,
        seed = Int(seed), total_ticks = 1000
    )
    for _ in 1:1000
        step!(model)
    end

    avg_deg = mean(degree(model.sharenet))
    endows = [a.harvest for a in allagents(model)|>collect]

    return (
        share_median_indegree = median(indegree(model.sharenet)),
        share_median_outdegree = median(outdegree(model.sharenet)),
        share_var_indegree = var(indegree(model.sharenet)),
        share_var_outdegree = var(outdegree(model.sharenet)),
        share_avg_clust = mean([a.share_clust for a in allagents(model)]),
        share_zero_indegree_count = count(==(0), indegree(model.sharenet)),
        pagerank_harvest_corr = avg_deg > 0 ? model.pagerank_harvest_corr = corspearman([a.pagerank for a in allagents(model)|>collect], endows) : 0.0,
        gini_indegree = gini(indegree(model.sharenet)),
    )
end

@everywhere get_summary_vector(sim) = Float64[
    sim.share_median_indegree,
    sim.share_median_outdegree,
    sim.share_var_indegree,
    sim.share_var_outdegree,
    sim.share_avg_clust,
    sim.share_zero_indegree_count,
    sim.pagerank_harvest_corr,
    sim.gini_indegree,
]

@everywhere function mahalanobis_distance(v, obs, inv_cov)
    sqrt((v .- obs)' * inv_cov * (v .- obs))
end

# latin–hypercube sampling from prior
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
    #seeds    = rand(rng, 1:10^7, N)
    # return a Vector of 7‐tuples
    [(l[i], C[i], γ[i], β[i], dens[i], ref[i], strength[i]) for i in 1:N]
end

@everywhere function log_prior(θ, rate)
    any(x->x<0, (θ[1],θ[2],θ[4])) && return -Inf
    any(x->x<0||x>1, (θ[3],θ[5],θ[6],θ[7])) && return -Inf
    log(rate[1]) + log(rate[2]) + log(rate[3]) - (rate[1]*θ[1] + rate[2]*θ[2] + rate[3]*θ[4])
end

@everywhere prior_pdf(θ, rate) = exp(log_prior(θ, rate))

# kernel functions
@everywhere tn_pos(μ,σ)=Truncated(Normal(μ,σ), 0, Inf)
@everywhere tn_01(μ,σ)=Truncated(Normal(μ,σ), 0, 1)

@everywhere function log_kernel_pdf(o, n, k::KernelBandwidth)
    lp = logpdf(tn_pos(o[1], k.σ_l), n[1])
    lp += logpdf(tn_pos(o[2], k.σ_C), n[2])
    lp += logpdf(tn_pos(o[4], k.σ_beta), n[4])
    lp += logpdf(tn_01(o[3], k.σ_γ), n[3])
    lp += logpdf(tn_01(o[5], k.σ_dens), n[5])
    lp += logpdf(tn_01(o[6], k.σ_ref),  n[6])
    lp += logpdf(tn_01(o[7], k.σ_strength), n[7])
    lp
end
@everywhere kernel_pdf(o,n,k) = exp(log_kernel_pdf(o,n,k))

@everywhere function kernel_sample(o, k::KernelBandwidth, rng)
    (
      rand(rng, tn_pos(o[1], k.σ_l)),
      rand(rng, tn_pos(o[2], k.σ_C)),
      rand(rng, tn_01(o[3], k.σ_γ)),
      rand(rng, tn_pos(o[4], k.σ_beta)),
      rand(rng, tn_01(o[5], k.σ_dens)),
      rand(rng, tn_01(o[6], k.σ_ref)),
      rand(rng, tn_01(o[7], k.σ_strength)),
      #rand(rng, 1:10^7)[1]
    )
end

# round 1: sample prior, estimate cov, then threshold
@everywhere function round1_fixedN_parallel!(
    obs_df::DataFrame,
    cols::Vector{Symbol},
    ε1::Float64,
    rate::Vector{<:Real},
    rng::AbstractRNG,
    N::Int
)
    # build observed summary
    obs = Float64[obs_df[1, c] for c in cols]

    # pilot LHS → one‐run sims to get covariance
    pilot_thetas = sample_from_prior_lhs(5*N, rate, rng)
    pilot_tasks = [(θ, rand(rng, UInt64)) for θ in pilot_thetas]
    pilot_summaries = @showprogress "[Round 1] pilot sims for cov…" pmap(pilot_tasks) do (θ, seed0)
        local_rng = Xoshiro(seed0)
        seed2 = rand(local_rng, 1:10^7)
        l, C, γ, β, dens, ref, strength = θ
        sim = run_model_expt((l, C, γ, β, dens, ref, strength, seed2))
        get_summary_vector(sim)
    end

    Mmat = hcat(pilot_summaries...)'
    @info "  [debug] M: typeof=$(typeof(Mmat)), size=$(size(Mmat))"
    Cmat = cov(Mmat)
    @info "  [debug] Cmat: typeof=$(typeof(Cmat)), size=$(size(Cmat))"
    invC = inv(Cmat + 1e-6I)
    @info "  [debug] invC: typeof=$(typeof(invC)), size=$(size(invC))"

    # accept‐N loop (batch‐parallel)
    acc_thetas = NTuple{7,Float64}[]
    acc_dists = Float64[]
    acc_hits = Int[]

    total_sims = 0
    while length(acc_thetas) < N
        need = N - length(acc_thetas)
        batch_size = max(ceil(Int, need * 1.2), need, 10)
        thetas = sample_from_prior_lhs(batch_size, rate, rng)
        tasks = [(θ, rand(rng, UInt64)) for θ in thetas]

        dmat = @showprogress "  evaluating $batch_size candidates…" pmap(tasks) do (θ, seed0)
            local_rng = Xoshiro(seed0)
            d = Vector{Float64}(undef, R_REPLICATES)
            l, C, γ, β, dens, ref, strength = θ
            for r in 1:R_REPLICATES
                seed2 = rand(local_rng, 1:10^7)
                sim = run_model_expt((l, C, γ, β, dens, ref, strength, seed2))
                d[r] = mahalanobis_distance(get_summary_vector(sim), obs, invC)
            end
            d
        end

        total_sims += R_REPLICATES * batch_size

        for (i, drow) in enumerate(dmat)
            length(acc_thetas) ≥ N && break
            md = median(drow)
            if md ≤ ε1
                push!(acc_thetas, thetas[i])
                push!(acc_dists,  md)
                push!(acc_hits,   count(x->x ≤ ε1, drow))
            end
        end
        @info "Round 1: accepted $(length(acc_thetas)) / $N so far"
    end

    # build weights from hit counts
    #w = Float64.(acc_hits)
    #w ./= sum(w)
    w = fill(1.0 / N, N)

    ess = 1/sum(w.^2)
    @info("Total sims = $total_sims")

    # return DataFrame of size N
    results = DataFrame(
      l = getindex.(acc_thetas, 1),
      C = getindex.(acc_thetas, 2),
      γ = getindex.(acc_thetas, 3),
      beta = getindex.(acc_thetas, 4),
      dens = getindex.(acc_thetas, 5),
      ref = getindex.(acc_thetas, 6),
      strength = getindex.(acc_thetas, 7),
      distance = acc_dists,
      hit_count = acc_hits,
      weight = w,
      ess = repeat([ess/N], length(w)),
      total_sims = repeat([total_sims], length(w))
    )

    # return BOTH df and invC
    return results, invC, total_sims
end

# subsequent rounds
@everywhere function round_t_fixedN_parallel!(
    prev_df::DataFrame,       # previous round particles
    prev_w::Vector{Float64},  # previous weights
    obs::Vector{Float64},
    invC::Matrix{Float64},
    rng::AbstractRNG,
    N::Int,                   # target accepted count
    α::Float64,               # quantile fraction
    rate::Vector{<:Real},
    kbw::KernelBandwidth,
    t::Int
)
    # fixed threshold from lagged quantile
    ε = quantile(prev_df.distance, α)
    @info "Round $t using ε = $ε"

    # prep for resampling
    cw = cumsum(prev_w)
    acc_thetas = NTuple{7,Float64}[]
    acc_dists  = Float64[]
    acc_hits   = Int[]

    # loop until we have N accepted
    n_acc = 0
    total_sims = 0

    while n_acc < N
      # how many more we need
      need = N - n_acc
      # choose a batch size somewhat larger than need/α to avoid repeated tiny batches
      batch_size = max(ceil(Int, need/α * 1.2), need, 10)

      # propose `batch_size` new thetas by resampling+kernel
      thetas = NTuple{7,Float64}[]
      for i in 1:batch_size
        #j = searchsortedfirst(cw, rand(rng))
        j = sample(rng, 1:length(prev_w), Weights(prev_w))
        θ_old = (
          prev_df.l[j], prev_df.C[j], prev_df.γ[j], prev_df.beta[j],
          prev_df.dens[j], prev_df.ref[j], prev_df.strength[j]
        )
        push!(thetas, kernel_sample(θ_old, kbw, rng))
      end

      # in parallel, simulate R replicates and compute distances
      tasks = [(θ, rand(rng, UInt64)) for θ in thetas]
      dmat = @showprogress "  evaluating $batch_size candidates…" pmap(tasks) do (θ, seed0)
        # each worker: one Xoshiro per candidate
        local_rng = Xoshiro(seed0)
        d = Vector{Float64}(undef, R_REPLICATES)
        for r in 1:R_REPLICATES
          seed2 = rand(local_rng, 1:10^7)
          l, C, γ, β, dens, ref, strength = θ
          sim = run_model_expt((l, C, γ, β, dens, ref, strength, seed2))
          d[r] = mahalanobis_distance(get_summary_vector(sim), obs, invC)
        end
        d
      end

      total_sims += R_REPLICATES * batch_size

      # filter accepted among this batch
      for (i, drow) in enumerate(dmat)
        if length(acc_thetas) >= N
          break
        end
        md = median(drow)
        if md ≤ ε
          push!(acc_thetas, thetas[i])
          push!(acc_dists, md)
          push!(acc_hits, count(x->x≤ε, drow))
          n_acc += 1
        end
      end

      @info "Round $t: accepted $(length(acc_thetas)) / $N so far"

    end

    # compute importance weights for the N accepted
    #old_nt = NamedTuple.(eachrow(prev_df))
    old_params = [(prev_df.l[i], prev_df.C[i], prev_df.γ[i], prev_df.beta[i],
               prev_df.dens[i], prev_df.ref[i], prev_df.strength[i])
              for i in 1:nrow(prev_df)]
    
    w = Float64[]
    for i in 1:N
      θ = acc_thetas[i]
      num = prior_pdf(θ, rate) #* acc_hits[i]
      #den = sum(prev_w[j]*kernel_pdf(old_nt[j], θ, kbw) for j in 1:length(prev_w))
      θ7 = (θ[1], θ[2], θ[3], θ[4], θ[5], θ[6], θ[7])
      den = sum(prev_w[j] * kernel_pdf(old_params[j], θ7, kbw) for j in 1:length(prev_w))
      push!(w, num/den)
    end
    w ./= sum(w)
    ess = 1/sum(w.^2)
    @info("Round $t ESS = $(ess)")
    @info("Total sims = $total_sims")

    # return DataFrame of size N
    results = DataFrame(
      l = getindex.(acc_thetas, 1),
      C = getindex.(acc_thetas, 2),
      γ = getindex.(acc_thetas, 3),
      beta = getindex.(acc_thetas, 4),
      dens = getindex.(acc_thetas, 5),
      ref = getindex.(acc_thetas, 6),
      strength = getindex.(acc_thetas, 7),
      distance = acc_dists,
      hit_count = acc_hits,
      weight = w,
      ess = repeat([ess/N], length(w)),
      total_sims = repeat([total_sims], length(w))
    )

    return results, total_sims
end


# driver
function abc_smc(obs_df, cols;
                ε1=3.0, nrounds=20, N=1000, α=0.3, 
                rate=[1.0, 0.1, 1.0], seed=80085,
                bandwidth_rule::Symbol = :gelman)

    rng = Xoshiro(seed)
    isfile(OUTFILE) && rm(OUTFILE)
    results = DataFrame[]

    #df1, invC = round1!(obs_df, cols, rate, rng, N, α)
    df1, invC, tsims = round1_fixedN_parallel!(
      obs_df, cols,
      ε1,     # your hand‐chosen starting threshold
      rate,   # e.g. 0.1
      rng,
      N       # desired number of particles
    )
    df1[!, :round] = fill(1, nrow(df1))
    push!(results, df1)
    CSV.write(OUTFILE, df1)
    @info "After Round 1: typeof(invC) = $(typeof(invC)), size = $(size(invC))"

    for t in 2:nrounds
        @info "=== Starting Round $t ==="
        prev = results[t-1]
        pw   = prev.weight
        # --- Adaptive bandwidth selection ---
        param_cols = [:l, :C, :γ, :beta, :dens, :ref, :strength]
        D  = length(param_cols)
        σs = [ std(prev[!, c]) for c in param_cols ]   # marginal std’s
        scalers = if bandwidth_rule == :silverman
        # Silverman multivariate rule-of-thumb
            h2 = (4.0/(D+2))^(2.0/(D+4)) * N^(-2.0/(D+4))
            sqrt(h2) .* σs
        elseif bandwidth_rule == :gelman
        # Gelman optimal MCMC scaling
            scale = 2.38 / sqrt(D)
            scale .* σs
        else
            error("Unknown bandwidth_rule: $bandwidth_rule")
        end
        kbw = KernelBandwidth(scalers...)

        df, tsims_t = round_t_fixedN_parallel!(
            prev, pw,
            Float64[obs_df[1, c] for c in cols],
            invC, rng, N, α, rate, kbw, t
        )

        tsims += tsims_t

        if isempty(df)
            @warn "Particle set depleted at round $t"
            break
        end

        @info "Grand total sims = $tsims"


        df[!, :round] = fill(t, nrow(df))
        push!(results, df)

        if mean(df[df.round .== t, :].ess) < 0.3 
            @warn "ESS below 0.3. Terminating...."
            break 
        end

        CSV.write(OUTFILE, df; append=true)

        if N / tsims_t < 0.005
            @warn "Acceptance rate below 0.5%. Terminating..."
            break
        end
    end

    return results
end

# let's go
col_list = [
  :share_median_indegree, 
  :share_median_outdegree,
  :share_var_indegree,
  :share_var_outdegree, 
  :share_avg_clust,
  :share_zero_indegree_count,
  :pagerank_harvest_corr,
  :gini_indegree,
]

obs_stats = CSV.read("../netdata/observed_stats.csv", DataFrame)
rounds = abc_smc(
    obs_stats, 
    col_list;
    ε1=3.0,
    nrounds=10, 
    N=2000, 
    α=0.2, 
    rate=[1.0, 0.1, 1.0], 
    #seed=6141,
    #seed=8008135,
    #seed=7457,
    seed=1993, 
    #seed=1337
    )

for (i, df) in enumerate(rounds)
    ess = 1 / sum(df.weight .^ 2)
    @info "Round $i: accepted=$(nrow(df)), ESS=$ess, dist=$(minimum(df.distance))--$(maximum(df.distance))"
end

CSV.write(OUTFILE, vcat(rounds...))
println("Done!")