### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 887a17bf-76fd-445b-a82d-ff5c5e45d3f9
begin
    using Pkg
    Pkg.activate("..")
    using Revise
	using SharingSignals
    using Agents, Graphs, Statistics, StatsBase, Distributions, Random
	using Plots, GraphRecipes, DataFrames, CSV, LaTeXStrings, KernelDensity, Colors, RecipesBase

	md"""
	## The emergence of sharing networks under indirect signaling
	#### Alejandro Pérez Velilla & Elspeth Ready
	"""
end

# ╔═╡ a5a9ecc8-5750-40fe-b433-e7f21d5b8fc8
begin
	kstar(u; c=6, beta=1.5, N=100) = clamp( floor( 1 + log(c/(u*N))/log( 1 - u*(1 - exp(-beta)) ) ), 0, Inf)
	
	kstar_plot = plot(
		0:0.01:1, kstar.(0:0.01:1), 
		lw=2, c="black", legend=true, legendtitle=L"\beta", label="", grid=false,
		xlabel=L"\mathrm{network\ density\ (u)}",
		ylabel=L"\mathrm{optimal\ number\ of\ sharing\ connections\ } (\hat{k})",
		xticks=(0.0:0.25:1.0, [L"%$a" for a in 0.0:0.25:1.0]),
		yticks=(0:5:25, [L"%$a" for a in 0:5:25])
	)
	scatter!(
		0:0.01:1, kstar.(0:0.01:1), 
		label=L"1.5", markercolor="white"
	)
	plot!(
		0:0.01:1, kstar.(0:0.01:1, beta=0.5), 
		lw=2, c="black", label=""
	)
	scatter!(
		0:0.01:1, kstar.(0:0.01:1, beta=0.5), 
		markercolor="gray", label=L"0.5"
	)
	plot!(
		0:0.01:1, kstar.(0:0.01:1, beta=0.25), 
		lw=2, c="black", label=""
	)
	scatter!(
		0:0.01:1, kstar.(0:0.01:1, beta=0.25), 
		markercolor="black", label=L"0.25"
	)

	savefig(kstar_plot, "../images/fig2_kstar.pdf")

	kstar_plot
end

# ╔═╡ da8d2640-88d5-406a-bf74-32575c229963
begin
	function lhs_sample_uniform(nsamples, rng)
	        step  = 1.0 / nsamples
	        pts   = step .* (collect(0:nsamples-1) .+ rand(rng, nsamples))
	        shuffle!(rng, pts)
	        return pts
	end
	
	lhs_sample(n,lo,hi,rng) = lo .+ (hi-lo).*lhs_sample_uniform(n,rng)
	
	lhs_sample_exponential(n,λ,rng) = -log.(lhs_sample_uniform(n,rng))./λ
	
	function sample_from_prior_lhs(N,rate,rng)
	        l = lhs_sample_exponential(N,rate[1],rng)
			C = lhs_sample_exponential(N,rate[2],rng)
	        γ = lhs_sample(N,0,1,rng)
			beta = lhs_sample_exponential(N,rate[3],rng) 
	        dens = lhs_sample(N,0,1,rng)
			ref = lhs_sample(N,0,1,rng)
			strength = lhs_sample(N,0,1,rng)
	        seeds = rand(rng,1:10^7,N)
	        [
				(l[i],C[i],γ[i],beta[i],dens[i],ref[i],strength[i],seeds[i]) 
				for i in 1:N
			]
	end

	round1_params = sample_from_prior_lhs(10000, [1.0, 0.1, 1.0], Xoshiro(4554))
	prior = DataFrame(
        l = [p[1] for p in round1_params],
        C = [p[2] for p in round1_params],
        γ = [p[3] for p in round1_params],
        beta = [p[4] for p in round1_params],
        dens = [p[5] for p in round1_params],
        ref = [p[6] for p in round1_params],
        strength = [p[7] for p in round1_params]
    )
	
	function plot_param_kde(
		posterior_df, prior_data, col;
	    title::String="", titlefontsize::Int=10,
		xlim::Tuple= (0,5), prior_color::String="black",
		prior_alpha::Float64=0.5, post_color="black", 
		xlab::AbstractString="", xlabelfontsize::Int=10, 
		rounds::Int64=1, pal::Symbol=:matter, xticks=false, yticks=false
	)
		posterior_data = posterior_df[!, col]
		
		c = cgrad(pal, rounds+1, categorical = true)
	    # Compute the KDEs for the posterior and prior samples
	    prior_kde = kde(prior_data)
	    
	    # Create the plot for the posterior KDE
	    p = plot(
			prior_kde.x, prior_kde.density,
			title=title, titlefontsize=titlefontsize,
			xlim=xlim, xlab=xlab, xlabelfontsize=xlabelfontsize,
			color=c[2], alpha=prior_alpha, lw=2, xticks=xticks, yticks=yticks
		)
		# Overlay the prior KDE on the same plot
		for i in 2:rounds
			post_kde = kde(
				posterior_df[posterior_df.round .== i,:][!, col]
			)
		    plot!(
				p, post_kde.x, post_kde.density,
		        color=c[i+1], lw=2, legend=false
			)
		end
	    return p
	end

	function plot_posterior_samples(
		df_post::DataFrame,
		df_prior::DataFrame,
		specs;
		ncols::Int = 2,
		legend::Bool = false,
		grid::Bool = false,
		size::Tuple{Int,Int} = (600,500),
		rounds::Int = 1,
		pal::Symbol = :matter,
	)

	    plots = [
	        plot_param_kde(
	            df_post,
	            df_prior[!, s.col],
				s.col;
	            xlab=s.xlab,
	            xlim=s.xlim,
	            titlefontsize=s.titlefontsize,
				xticks=s.xticks,
				yticks=s.yticks,
				rounds=rounds,
				pal=pal
	        )
	        for s in specs
	    ]
	
	    n = length(plots)
	    nrows = ceil(Int, n / ncols)
	    return plot(plots...;
	                layout=(nrows, ncols),
	                legend=legend,
	                grid=grid,
	                size=size,
					dpi=300
				   )
	end

	gaussian_kernel(d, scale) = @. exp(-0.5 * (d/scale)^2)

	ess(w) = 1 / sum(abs2, w)

	function normalize!(w::AbstractVector{<:Real})
	    s = sum(w)
	    s <= 0 && error("Sum of weights <= 0; cannot normalize.")
	    w ./= s
	    return w
	end
	
	function aggregate_posts(dfs)
		posts = []
		for d in dfs
			post = CSV.read(d, DataFrame)
			#post = post[post.round .== maximum(post.round), :]
			push!(posts, post)
		end
		pooled = vcat(posts...)
        return pooled
	end

	function normalize_pooled_post(pooled; epsilon_star=1.0, q_star = 0.8, temp = 1.0)
		# absolute-distance kernel
        dists = pooled.distance
        eps_used = isnothing(epsilon_star) ? quantile(dists, q_star) : epsilon_star
        scale = eps_used * temp
        K = gaussian_kernel(dists, scale)

        pooled.weight_old = copy(pooled.weight)
        pooled[!, :weight] = pooled[!, :weight] .* K
        normalize!(pooled[!, :weight])
		return pooled
	end

	col_list = [
	  :share_avg_degree, 
	  :share_median_indegree, 
	  :share_median_outdegree,
	  :share_avg_clust,  
	  :share_coreness_avg,     
	  :share_var_indegree,
	  :share_var_outdegree, 
	  :share_zero_indegree_count,
	  :pagerank_harvest_corr,
	  :gini_indegree
	]

	specs = [
	    (
			col=:c, 
			xlab=L"\mathrm{connection\ cost\ scale\ } (c)", 
			xlim=(0,60), titlefontsize=12, 
			xticks=(0:10:60, [L"%$a" for a in 0:10:60]), yticks=false
		),
	    (
			col=:l, xlab=L"\mathrm{connection\ cost\ elasticity\ } (l)", 
			xlim=(0,10), titlefontsize=12, 
			xticks=(0:2:10, [L"%$a" for a in 0:2:10]), yticks=false
		),
	    (
			col=:gamma, xlab=L"\mathrm{harvest\ advantage\ } (γ)", 
			xlim=(0,1), titlefontsize=12, 
			xticks=(0:0.2:1.0, [L"%$a" for a in 0:0.2:1.0]), yticks=false
		),
	    (
			col=:beta, xlab=L"\mathrm{social\ reward\ } (β)", 
			xlim=(0,3), titlefontsize=12, 
			xticks=(0:0.5:3, [L"%$a" for a in 0:0.5:3]), yticks=false
		),
	    (
			col=:dens, xlab=L"\mathrm{comm.\ network\ density\ } (u)", 
			xlim=(0,1),  titlefontsize=12, 
			xticks=(0:0.2:1.0, [L"%$a" for a in 0:0.2:1.0]), yticks=false
		),
	    (
			col=:ref, xlab=L"\mathrm{central\ harvest\ } (\zeta)", 
			xlim=(0,1), titlefontsize=12, 
			xticks=(0:0.2:1.0, [L"%$a" for a in 0:0.2:1.0]), yticks=false
		),
	    (
			col=:strength, xlab=L"\mathrm{correlation\ strength\ } (S)", 
			xlim=(0,1), titlefontsize=12, 
			xticks=(0:0.2:1.0, [L"%$a" for a in 0:0.2:1.0]), yticks=false
		),
	]
	
	# split them into “preference” vs “comm”:
	pref_specs = specs[1:4]
	comm_specs = specs[5:7]

	seed1 = "../data/mad_fit1.csv"
	seed2 = "../data/mad_fit2.csv"
	seed3 = "../data/mad_fit3.csv"
	seed4 = "../data/mad_fit4.csv"
	seed5 = "../data/mad_fit5.csv"
	
	post_unf = aggregate_posts([seed1, seed2, seed3, seed4, seed5])
	#CSV.write("../data/mad_fit_pooled.csv", post_unf)
	
	post_unf.c .= post_unf.C
	post_unf.gamma .= post_unf.γ
	post = post_unf[post_unf.round .== maximum(post_unf.round),:]
	post.round .= repeat([2], nrow(post))
	post = normalize_pooled_post(post)
	rounds = maximum(post.round)
	prior = CSV.read("../data/prior_pred.csv", DataFrame)
	
	#each a single‐column plot
	p_pref = plot_posterior_samples(post, prior, pref_specs; ncols=1, rounds=rounds, pal=:binary)
	p_comm = plot_posterior_samples(post, prior, comm_specs; ncols=1, rounds=rounds, pal=:binary)
	
	#tile them side by side
	postplot = plot(p_pref, p_comm; layout=(1,2), size=(800,600))

	savefig(postplot, "../images/fig3_postplot.pdf")

	postplot
end

# ╔═╡ 13a5a8bc-e920-46c8-b24e-0fee93f4cf8b
ess_post = ess( post.weight ) / nrow(post)

# ╔═╡ 2c6e7505-8036-4a84-bad8-19d861034a53
md"""
#### Tables
"""

# ╔═╡ 28f2262c-56e8-429c-81e9-ecdc432cd691
begin
	# map unicode → LaTeX for Greek letters
	const _GREEK = Dict(
	    "α" => "\\alpha",  "β" => "\\beta",   "γ" => "\\gamma",
	    "Δ" => "\\Delta",  "μ" => "\\mu",     "σ" => "\\sigma",
		"ζ" => "\\zeta"
	)
	
	# always wrap in math‐mode, but substitute macros where known
	function latex_label(sym::Symbol)
	    s = string(sym)
	    mac = get(_GREEK, s, s)
	    return "\$" * mac * "\$"
	end
	
	"""
	    df_to_latex_stats(df, cols; digits=3)
	
	For each column in `cols` of `df`, compute:
	  • Mean  
	  • Median  
	  • Std  
	  • 75% CI  
	  • 87% CI  
	  • 99% CI  
	
	rounding all to `digits`, and return a booktabs‐style LaTeX `tabular` as a String.
	"""
	function df_to_latex_stats(df::DataFrame, cols::Vector{Symbol}; digits::Int=3)
	    cis = [0.75, 0.87, 0.95]
	    ci_names = ["$(Int(round(p*100)))\\% CI" for p in cis]
	
	    rows = Vector{Vector{String}}()
	    for c in cols
	        x = skipmissing(df[!, c])
	        μ = round(mean(x),   digits=digits)
	        med = round(median(x), digits=digits)
	        σ = round(std(x),    digits=digits)
	
	        ci_strs = String[]
	        for p in cis
	            lo = round(quantile(x, (1-p)/2), digits=digits)
	            hi = round(quantile(x, (1+p)/2), digits=digits)
	            push!(ci_strs, "\$[$lo,\\;$hi]\$")
	        end
	
	        push!(rows,
	            [ latex_label(c),
	              string(μ),
	              string(med),
	              string(σ),
	              ci_strs...  
	            ]
	        )
	    end
	
	    ncols = 4 + length(cis)
	    colspec = "l" * join(fill("c", ncols-1), "")
	
	    header = "\\begin{tabular}{" * colspec * "}\n" *
	             "\\toprule\n" *
	             "Parameter & Mean & Median & Std & " *
	             join(ci_names, " & ") * " \\\\\n" *
	             "\\midrule\n"
	
	    # insert \addlinespace after each data row
	    body = ""
	    for row in rows
	        body *= join(row, " & ") * " \\\\\n"
	        body *= "\\addlinespace\n"
	    end
	
	    footer = "\\bottomrule\n\\end{tabular}"
	
	    return header * body * footer
	end


	#post = post_smcr[post_smcr.round .== rounds, :]
	post.β = post.beta
	post.c = post.C
	post.u = post.dens
	post.ζ = post.ref
	post.S = post.strength
	
	println(df_to_latex_stats(post[post.round .== maximum(post.round),:], [:c, :l, :γ, :β, :u, :ζ, :S]))
end

# ╔═╡ 85943361-ab1e-4266-9c59-c4184b49948f
md"""
### Supplementary figures
"""

# ╔═╡ 8db4d3f1-5060-46a6-a9ff-4ca2941b5247
begin
	data_net, harvest, data_net_unfiltered = datanet()

	#netdatas = DataFrame(harvest=harvest, indegree=indegree(data_net), outdegree=outdegree(data_net), clustering=local_clustering_coefficient(data_net), pagerank=pagerank(data_net))

	data_stats = extract_netdata_full(data_net, log.(1 .+ harvest./2500))
	#CSV.write("../data/observed_stats.csv", data_stats)

	#data_stats_unf = extract_netdata_full(data_net_unfiltered, log.(1 .+ harvest./2500))
	#CSV.write("../data/observed_stats_unf.csv", data_stats_unf)
	
	default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
	
	harvest_plot = plot(
		histogram(
			harvest./2500, 
			xlab=L"\mathrm{harvest\ (kCal/D)}",
			bins=50,
			legend=false, 
			grid=false, 
			xlabelfontsize=12,
			xticks=(0:1000:4000, [L"%$a" for a in 0:1000:4000]),
			yticks=(0:10:60, [L"%$a" for a in 0:10:60])
		),
		histogram(
			log.(1 .+ harvest./2500), 
			bins=50, 
			xlab=L"\mathrm{harvest\ (kCal/D,\ } \log(1+x))", 
			legend=false, 
			grid=false, 
			xlabelfontsize=12,
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
			yticks=(0:10:30, [L"%$a" for a in 0:10:30])
		),
		size=(500, 300), dpi=300
	)

	#savefig(harvest_plot, "../images/sup1_harvest.pdf")
end

# ╔═╡ e22c7a83-bff9-4766-bcf4-071077150b13
med_harvest = median(harvest ./ 2500)

# ╔═╡ cc9089be-7dc3-400a-91dc-1f5c442ee41d
ecdf(harvest ./ 2500)(0)

# ╔═╡ 988cd7d4-b6f4-4579-bdca-45ee88854f38
begin
	post_pred = CSV.read("../data/post_pred_mad.csv", DataFrame)
	replace!(prior.reciprocity_dyad, missing => 0.0)
	replace!(prior.reciprocity_dyad, NaN => 0.0)
	replace!(prior.reciprocity_edge, missing => 0.0)
	replace!(prior.reciprocity_edge, NaN => 0.0)

	col_list1 = [
	    :share_avg_degree,
	    :share_median_outdegree,
		:share_median_indegree,
		:share_coreness_avg,
		#:share_zero_indegree_count,
	]

	col_list2 = [
		:share_avg_clust,
		:gini_indegree,
		:reciprocity_dyad,
		:reciprocity_edge
	]

	col_list3 = [
	 	:share_var_outdegree,
		:share_var_indegree
	]

	col_list4 = [
		:pagerank_harvest_corr,
	]
	
	function plot_prior_post(
		prior_df::DataFrame, post_df::DataFrame, 
		data_df::DataFrame, cols::Vector{Symbol};
	    q_low=0.05, q_high=0.95, offset=0.2, ylabel=L"\mathrm{value}", 
	    xlim = (0.5, length(cols)+0.5), 
		ylim = (-0.1, 90),
	    legend = true, 
		xnames = (), 
		ynames = (
			[0, 10, 20, 30, 40, 50, 60, 70, 80], 
			[L"%$a" for a in [0, 10, 20, 30, 40, 50, 60, 70, 80]]
		)
	)
	    # Prepare vectors to store summary statistics
	    prior_means = Float64[]
	    prior_lower_errors = Float64[]
	    prior_upper_errors = Float64[]
	
	    post_means = Float64[]
	    post_lower_errors = Float64[]
	    post_upper_errors = Float64[]
	
	    data_means = Float64[]
	
	    # Compute summary statistics for each specified column
	    for col in cols
	        # Prior statistics (using median as the central measure)
	        data_prior = collect(skipmissing(prior_df[!, col]))
	        m_prior = median(data_prior)
	        lower_prior = quantile(data_prior, q_low)
	        upper_prior = quantile(data_prior, q_high)
	        push!(prior_means, m_prior)
	        push!(prior_lower_errors, m_prior - lower_prior)
	        push!(prior_upper_errors, upper_prior - m_prior)
	
	        # Post statistics (using median as the central measure)
	        data_post = collect(skipmissing(post_df[!, col]))
	        m_post = median(data_post)
	        lower_post = quantile(data_post, q_low)
	        upper_post = quantile(data_post, q_high)
	        push!(post_means, m_post)
	        push!(post_lower_errors, m_post - lower_post)
	        push!(post_upper_errors, upper_post - m_post)
	
	        # Data statistics (using the mean)
	        data_data = collect(skipmissing(data_df[!, col]))
	        m_data = mean(data_data)
	        push!(data_means, m_data)
	    end
	
	    ncols = length(cols)
	    # Define x positions: prior is shifted left, post is shifted right, and data is centered.
	    x_prior = [i - offset for i in 1:ncols]
	    x_post  = [i + offset for i in 1:ncols]
	    x_data  = 1:ncols
	    xticks = isempty(xnames) ? (1:ncols, string.(cols)) : xnames
	    yticks = ynames
	
	    # Create the scatter plot with error bars for the prior data
	    p = scatter(
	        x_prior, prior_means,
	        yerror = (prior_lower_errors, prior_upper_errors),
	        marker = (:circle, 4),
	        label = L"\mathrm{prior}",
	        ylabel = ylabel,
	        xticks = xticks,
	        yticks = yticks,
	        xrotation = 90,
	        xlim = xlim,
	        ylim = ylim,
	        legend = legend,
	        color = :black,
	        grid = false
	    )
	
	    # Overlay the post data with error bars
	    scatter!(
	        p, x_post, post_means,
	        yerror = (post_lower_errors, post_upper_errors),
	        marker = (:diamond, 4),
	        label = L"\mathrm{post}",
	        color = "black",
			fillcolor = :black
	    )
	
	    # Overlay the data mean (plotted without error bars) at the center
	    scatter!(
	        p, x_data, data_means,
	        marker = (:xcross, 4),
	        label = L"\mathrm{data}",
	        markercolor = :black
	    )
	
	    return p
	end

	function plot_post_pred(prior, post_pred, data_stats)
		return plot(
			plot_prior_post(
				prior, 
				post_pred, 
				data_stats,
				col_list1, 
				xnames = ([1,2,3,4],[L"\mathrm{MEAN.\ DEG}", L"\mathrm{MED.\ IN}", L"\mathrm{MED.\ OUT}", L"\mathrm{CORE}"]),
				ynames = (0:20:100, [L"%$a" for a in 0:20:100]),
				ylim = (-0.5, 100),
				ylabel = "",
				legend = :topleft
			),
			plot_prior_post(
				prior, 
				post_pred,
				data_stats,
				col_list2, 
				xnames = ([1, 2, 3, 4], [L"\mathrm{CLUST}", L"\mathrm{GINI}", L"\mathrm{REC. DYAD}", L"\mathrm{REC. EDGE}"]),
				ynames = (0.0:0.2:0.8, [L"%$a" for a in 0.0:0.2:0.8]),
				ylabel = "", 
				offset = 0.2, 
				xlim = (0.5, 4.5), 
				ylim = (-0.01, 0.8),
				legend = false
			),
			plot_prior_post(
				prior, 
				post_pred,
				data_stats,
				col_list3, 
				xnames = ([1, 2],[L"\mathrm{VAR.\ IN}", L"\mathrm{VAR.\ OUT}"]),
				ynames = (0:100:600, [L"%$a" for a in 0:100:600]),
				ylabel="", 
				offset=0.2, 
				xlim=(0.5, 2.5), 
				ylim = (-0.5, 600),
				legend=false
			),
			plot_prior_post(
				prior, 
				post_pred,
				data_stats,
				col_list4, 
				xnames = ([1],[L"\mathrm{PGRNK}"]),
				ynames = (-1.0:0.2:1.0, [L"%$a" for a in -1.0:0.2:1.0]),
				ylabel="", 
				offset=0.4, 
				xlim=(0.0, 2.5), 
				ylim = (-1, 1),
				legend=false
			),
			layout = @layout([a{0.4w} b{0.35w} c{0.15w} d{0.1w}]), bottom_margins=8Plots.mm, left_margins=1Plots.mm,size=(800, 400), dpi=300
		)
	end

	post_pred_plot = plot_post_pred(prior, post_pred, data_stats)

	savefig(post_pred_plot, "../images/fig3_postpred.pdf")

	post_pred_plot
end

# ╔═╡ d8a0ccfc-c28f-446c-b61d-e0ca8933dec8
begin
	function print_prob_table(prior, post; digits=3)
	    # compute the three quantities
	    p1 = mean(   prior.com_endow_corr .< 0)
	    p2 = mean(   prior.pagerank_harvest_corr .< 0)
	    p3 = mean(( prior.com_endow_corr .< 0) .& (prior.pagerank_harvest_corr .< 0))
	
	    q1 = mean(   post.com_endow_corr .< 0)
	    q2 = mean(   post.pagerank_harvest_corr .< 0)
	    q3 = mean(( post.com_endow_corr .< 0) .& (post.pagerank_harvest_corr .< 0))
	
	    # round
	    p1,p2,p3 = round.( (p1,p2,p3), digits=digits )
	    q1,q2,q3 = round.( (q1,q2,q3), digits=digits )
	
	    println("\\begin{table}[htbp]")
	    println("  \\centering")
	    println("  \\begin{tabular}{lcc}")
	    println("    \\toprule")
	    println("    Metric & Prior & Posterior \\\\")
	    println("    \\midrule")
	    println("    \$\\Pr(\\mathrm{comm.\\ endow.\\ corr.\\ } < \\ 0)\$ & $p1 & $q1 \\\\")
		println("    \\addlinespace")
	    println("    \$\\Pr(\\mathrm{share\\ PageRank\\ endow.\\ corr.\\ } < \\ 0)\$ & $p2 & $q2 \\\\")
		println("    \\addlinespace")
	    println("    \$\\Pr(\\mathrm{both\\ } <\\ 0)\$ & $p3 & $q3 \\\\")
	    println("    \\bottomrule")
	    println("  \\end{tabular}")
	    println("  \\caption{Prior vs. posterior probabilities}")
	    println("  \\label{tab:prior_post_probs}")
	    println("\\end{table}")
	end
	
	print_prob_table(prior, post_pred)
end

# ╔═╡ b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
begin
	post_end = post[post.round .== maximum(post.round),:]
	#post_end = varpost[varpost.round .== maximum(varpost.round),:]

	param_cols = [:l,:C,:γ,:beta,:dens,:ref,:strength]
	wv = collect(post_end.weight); wv ./= sum(wv)
	μ = [sum(post_end[!,c] .* wv) for c in param_cols]
	θs = [(post_end.l[i], post_end.C[i], post_end.γ[i], post_end.beta[i],
	       post_end.dens[i], post_end.ref[i], post_end.strength[i]) for i in 1:nrow(post_end)]
	i_medoid = argmin([sum((θs[i] .- μ).^2) for i in eachindex(θs)])
	θ = θs[i_medoid]

	
	
	model = initialize_sharing_signals_ywb(
	    C = θ[2], l = θ[1], γ = θ[3], beta = θ[4],
	    dens = θ[5], ref = θ[6], strength = θ[7],
	    total_ticks = 1000, seed = 745745
	)

	run!(
		model, 1001, 
		mdata = [
			:l,
			:h,
			:γ,
			:beta,
			:dens,
			:strength,
			:ref,
	        :share_density,
	        :com_density,
	        :share_avg_degree,
	        :com_avg_degree,
	        :share_avg_clust,
	        :com_avg_clust,
	        :share_giant_component,
	        :com_giant_component,
	        :share_diameter,
	        :com_diameter,
	        :share_avg_pathlength,
	        :com_avg_pathlength,
	        :com_degree_assort,
	        :share_degree_assort,
	        :com_wealth_assort,
	        :share_wealth_assort,
	        :com_endow_deg_corr,
	        :share_endow_indeg_corr,
	        :com_zero_deg_count,
	        :share_zero_outdegree_count,
	        :share_zero_indegree_count,
	        :com_endow_deg_corr_nonzero,
	        :share_endow_deg_corr_nonzero,
	        :com_coreness_avg,
	        :com_coreness_max,
	        :com_coreness_endow_corr,
	        :share_coreness_avg,
	        :share_coreness_max,
	        :share_coreness_endow_corr
    ],
		when_model = 1000
	)

	
	function net_hist(model, data_net; model_inbins=45, model_outbins=25)
		endhist = plot(
			plot(
				histogram(
					[a.indeg for a in allagents(model)|>collect], ylab=L"\mathrm{model}", 
					legend=false, 
					xlim=(-0.1, 30),
					ylim=(0,25),
					color=palette(:tokyo10)[3], 
					bins=model_inbins,
					grid=false,
					xticks=(0:5:25, [L"%$a" for a in 0:5:25]),
					yticks=(0:5:25, [L"%$a" for a in 0:5:25]),
				),
				histogram(
					indegree(data_net), 
					xlab=L"\mathrm{indegree}", 
					ylab=L"\mathrm{data}",
					xlim=(-0.1, 30),
					ylim=(0,25),
					legend=false, bins=25, 
					color=palette(:tokyo10)[7],
					grid=false,
					xticks=(0:5:25, [L"%$a" for a in 0:5:25]),
					yticks=(0:5:25, [L"%$a" for a in 0:5:25]),
				),
				layout=(2,1), link=:all
			),
			plot(
				histogram(
					[a.outdeg for a in allagents(model)|>collect], #xlab=L"\mathrm{outdegree}", 
					legend=false, xlim=(0,35), 
					ylim=(0, 40), bins=model_outbins, 
					color=palette(:tokyo10)[3],
					grid=false,
					xticks=(0:5:35, [L"%$a" for a in 0:5:35]),
					yticks=(0:5:40, [L"%$a" for a in 0:5:40]),
				),
				histogram(
					outdegree(data_net), 
					xlab=L"\mathrm{outdegree}", 
					legend=false, bins=30, 
					ylim=(0,40), 
					xlim=(0,35),
					color=palette(:tokyo10)[7],
					grid=false,
					xticks=(0:5:35, [L"%$a" for a in 0:5:35]),
					yticks=(0:5:40, [L"%$a" for a in 0:5:40]),
				),
				layout=(2,1), link=:all
			),
			layout=(1,2)
		) 
		
		end_out = scatter(
			[log(1 + a.harvest) for a in allagents(model)|>collect],
			[a.outdeg for a in allagents(model)|>collect],
			legend=false,
			#xlab=L"\mathrm{log\ endowment}",
			xticks=false,
			yticks=(0:5:25, [L"%$a" for a in 0:5:25]),
			ylab=L"\mathrm{outdegree}",
			title=L"\mathrm{model}",
			alpha=0.25,
			color=palette(:tokyo10)[3],
			grid=false,
			ylim=(-1,25),
			xlim=(-1,8)
		)

		end_in = scatter(
			[log(1 + a.harvest) for a in allagents(model)|>collect],
			[a.indeg for a in allagents(model)|>collect],
			#xlab=L"\mathrm{log\ harvest}",
			ylab=L"\mathrm{indegree}",
			yticks=(0:5:25, [L"%$a" for a in 0:5:25]),
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
			legend=false,
			alpha=0.25,
			color=palette(:tokyo10)[3],
			grid=false,
			ylim=(-1,25),
			xlim=(-1,8)
		)
	
		data_out = scatter(
			log.(1 .+ (harvest./2500)),
			outdegree(data_net),
			legend=false,
			#xlab=L"\mathrm{log\ endowment}",
			#ylab=L"\mathrm{outdegree}",
			title=L"\mathrm{data}",
			yticks=false,
			xticks=false,
			alpha=0.25,
			color=palette(:tokyo10)[7],
			grid=false,
			ylim=(-1,25),
			xlim=(-1,8)
		)
	
		data_in = scatter(
			log.(1 .+ (harvest./2500)),
			indegree(data_net),
			legend=false,
			#xlab=L"\mathrm{log\ harvest}",
			yticks=false,
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
			#ylab=L"\mathrm{indegree}",
			alpha=0.25,
			color=palette(:tokyo10)[7],
			grid=false,
			ylim=(-1,25),
			xlim=(-1,8)
		)
	
		sim_clust = scatter(
			[log(1 + a.harvest) for a in allagents(model)|>collect],
			[a.share_clust for a in allagents(model)|>collect],
			xlab=L"\log(1+\mathrm{harvest})",
			ylab=L"\mathrm{clustering}",
			legend=false,
			color=palette(:tokyo10)[3],
			alpha=0.5,
			xlim=(-1,8),
			ylim=(-0.01, 0.27),
			grid=false,
			yticks=(0.0:0.05:0.25, [L"%$a" for a in 0.0:0.05:0.25]),
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
		)
	
		data_clust = scatter(
			log.(1 .+ harvest./2500),
			local_clustering_coefficient(data_net),
			legend=false,
			xlab=L"\log(1+\mathrm{harvest)}",
			color=palette(:tokyo10)[7],
			alpha=0.5,
			xlim=(-1,8),
			ylim=(-0.01, 0.27),
			yticks=false,
			grid=false,
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
		)
	
		p2 = plot(
			end_out,
			end_in,
			sim_clust,
			layout=(3,1)
		)
	
		p3 = plot(
			data_out,
			data_in,
			data_clust,
			layout=(3,1)
		)
	
		p4 = plot(p2, p3, layout=(1,2))
	
		plot(
			plot(endhist, p4, layout=(1,2), size=(800, 400), margins=1.5Plots.mm, dpi=300),
			layout=(2,1)
		)
	end

	nethist_plot = net_hist(model, data_net)

	savefig(nethist_plot, "../images/fig4_nethist.pdf")

	nethist_plot
end

# ╔═╡ d4bb1701-efa7-4bb0-ab69-73e0b7b8f5ac
modal_zeta = θs[i_medoid][6]

# ╔═╡ 56af3b89-c9aa-405e-9deb-cabb77e69993
central_harvest = exp( (θs[i_medoid][6] .* maximum(log.(1 .+ harvest ./ 2500))) ) - 1

# ╔═╡ 88f07c6d-ad7f-4d9f-845b-3ef91cc59f28
begin
	function plot_coreness(g, h; title="")
		
		cn = core_number(g)
		
		# ── 1) Compute coreness & group into shells ────────────────
		groups = Dict{Int, Vector{Int}}()
		for v in vertices(g)
		    push!( get!(groups, cn[v], Int[]), v )
		end
		shell_keys = sort(collect(keys(groups)), rev=true)
		shells    = [ groups[k] for k in shell_keys ]
		
		# ── 2) Manually build concentric coordinates ───────────────
		n  = nv(g)
		coords = zeros(2, n)
		ns = length(shells)
		
		for (i, shell) in enumerate(shells)
		    # spread shells from r=1/(ns+1) up to ns/(ns+1)
		    r = i/(ns+1)
		
		    m = length(shell)
		    for (j, v) in enumerate(shell)
		        θ = 2π*(j-1)/m
		        coords[1, v] = r*cos(θ)
		        coords[2, v] = r*sin(θ)
		    end
		end
		
		# 3) build a discrete Viridis palette keyed to core-values
		minc, maxc = minimum(values(cn)), maximum(values(cn))
		nlevels = maxc - minc + 1
		pal = palette(:coffee, nlevels)   # Array of nlevels Colors
		
		# now map each node’s core to a color
		node_colors = [ pal[ cn[v] - minc + 1 ] for v in vertices(g) ]
		
		# unpack & plot
		x, y = coords[1, :], coords[2, :]
		graphplot(
		    g;
		    x = x,
		    y = y,
			curves = false,
			arrow = 0.1,
			shorten = 0.1,
			arrow_scale = 0.01,
			node_weights=h,
		    nodecolor = node_colors,
		    node_size = 0.3,
		    node_labels = vertices(g),
		    linecolor = :gray,
			title = title,
			dpi = 300,
			size = (800, 400),    # width×height in pixels
		)
	
	end
	
	function plot_core(g1, g2; g3=true, maxcore=9)
		coreness1 = []
		count = 0
		for i in 0:maxcore
			count += length(k_shell(g1, i))
			push!(coreness1, count)
		end
		count = 0
		coreness2 = []
		for i in 0:maxcore
			count += length(k_shell(g2, i))
			push!(coreness2, count)
		end
		if g3
			rand_graph = erdos_renyi(110, mean(indegree(g1))/110)
			count = 0
			coreness3 = []
			for i in 0:maxcore
				count += length(k_shell(rand_graph, i))
				push!(coreness3, count)
			end
		end
			
		data_coreplot = plot(
			0:maxcore, coreness1[1:end], 
			color="black", ls=:dot,
			lw=2, label=L"\mathrm{data}", 
			xlabel=L"\mathrm{core\ number}", ylabel=L"\mathrm{cumulative\ count}", 
			size=(500,500), 
			legend=:topleft, 
			xticks=(
				0:2:maxcore,
				[L"%$a" for a in 0:2:10]
			), 
			yticks=(
				0:20:100, 
				[L"%$a" for a in 0:20:100]
			)
		)
		plot!(
			0:maxcore,
			coreness2[1:end], 
			color="black", alpha=0.5,
			lw=2, label=L"\mathrm{model}"
		)
		if g3 plot!(0:maxcore, coreness3[1:end], color="black", ls=:dash, lw=2, label=L"\mathrm{random}") end
		scatter!(0:maxcore, coreness1[1:end], lw=2, label="")
		scatter!(0:maxcore, coreness2[1:end], lw=2, label="")
		if g3 scatter!(0:maxcore, coreness3[1:end], lw=2, label="") end
	end

	netplots = plot(
		plot_coreness(model.sharenet, (1000 .+ harvest).^1.1, title=L"\mathrm{model}",),
		plot_coreness(data_net, (1000 .+ harvest).^1.1, title=L"\mathrm{data}"),
		layout=(1,2), dpi=300
	)
	
	coreplot = plot(
		netplots,
		plot_core(data_net, model.sharenet),
		layout=(2,1),
		size=(600,600)
	)

	savefig(coreplot, "../images/fig5_coreplot.pdf")

	coreplot
end

# ╔═╡ a657563c-fa1e-455d-ac70-3fabdd52281c
begin
	function plot_corr(dat; title=L"\mathrm{A.\ prior}")
	    x, y = dat.com_endow_corr, dat.pagerank_harvest_corr
	
	    # 1) build the three real panels + one blank
	    hx = histogram(
			x;
			orientation = :v,
			framestyle = :none,
			legend = false,
			bins = 90,
			#title = title,
			titlelocation = :left,
			margin = 0Plots.mm,
			xlims = (-1,1),
	    )
	
	    s = scatter(
			x, y;
			framestyle = :box,
			legend = false,
			ylims = (-1, 1),
			xlims = (-1, 1),
			alpha = 0.05,
			#title = title,
			markerstrokewidth = 0,
			markercolor = :black,
			xlab = L"\mathrm{comm.\ deg\!-\!harvest\ corr}",
			ylab = L"\mathrm{share.\ PageRank\!-\!harvest\ corr}",
			xlabelfontsize = 14,
			ylabelfontsize = 14,
			xticks = (-1:0.5:1, [L"%$a" for a in -1:0.5:1]),
			yticks = (-1:0.5:1, [L"%$a" for a in -1:0.5:1]),
			margin = 0Plots.mm,
			dpi = 300,
	    )
	    vline!(s, [0], color=:black, alpha=1.0)
	    hline!(s, [0], color=:black, alpha=1.0)
		hline!([data_stats.pagerank_harvest_corr], lw=2, color=:grey, ls=:dot)
	
	    hy = histogram(
	      y;
	      orientation = :h,
	      framestyle = :none,
	      legend = false,
	      margin = 0Plots.mm,
	      ylims = (-1,1),
	    )
	
	    blank = plot(
	      [];
	      framestyle = :none,
	      xaxis = false,
	      yaxis = false,
	      margin = 0Plots.mm,
	    )
	
	    # 2) carve out a 2×2 grid: 
	    #    heights =  20% for top row,  80% for bottom
	    #    widths  =  80% for left col, 20% for right
	    l = RecipesBase.grid(2, 2; heights=(0.1,0.9), widths=(0.9,0.1))
	
	    # 3) place the panels in row-major order:
	    #    (1,1)=hx  (1,2)=blank
	    #    (2,1)=s   (2,2)=hy
	    plot(
	      	hx,  blank,
	      	s,   hy;
	      	layout = l,
	      	size = (500, 500),
	      	bottom_margin = 11Plots.px,
	      	left_margin = 7Plots.mm,
			plot_title = title,
			plot_titlelocation=:left
	    )
	end

	layout = @layout [a{0.475w} _{0.05w} b{0.475w}]

	postpred_nbt = plot(
		plot_corr(prior),
		vline([0], yaxis=false, xaxis=false, lw=3, color=:black),
		plot_corr(post_pred, title=L"\mathrm{B.\ post}"),
		layout=layout,
		size=(1200,600), dpi=300
	)

	savefig(postpred_nbt, "../images/fig6_nbt.pdf")

	postpred_nbt
end

# ╔═╡ 88fefd48-f055-4a0e-b56f-1f38d8eae9bc
begin
	left_scatter = scatter(
		log.( 1 .+ [a.harvest for a in allagents(model)|>collect] ),
		[length(a.comneighbors) for a in allagents(model)|>collect],
		ylab=L"\mathrm{comm.\ degree}",
		ylabelfontsize=17,
		xticks=(0:8, [L"%$a" for a in 0:8]),
		yticks=(0:10:50, [L"%$a" for a in 0:10:50]),
		alpha=0.5, ylim=(0,50),
		markercolor=:black,
		xlim=(-0.5, 9)
	)
	annotate!([10], [-5.5], text(L"\log (1 + \mathrm{harvest})", 17))

	right_scatter = plot(
		scatter(
			log.( 1 .+ [a.harvest for a in allagents(model)|>collect] ),
			[a.pagerank for a in allagents(model)|>collect],
			#[a.indeg for a in allagents(model)|>collect],
			ylab=L"\mathrm{share\ PageRank\ (model)}",
			ylabelfontsize=8,
			xticks=(0:8, [L"%$a" for a in 0:8]),
			alpha=0.5,
			markercolor=:black,
			yticks=(0:0.005:0.03, [L"%$a" for a in 0:0.005:0.03]),
			ylim=(0, 0.035),
			xlim=(-0.5, 9)
		),
		scatter(
			log.(1 .+ harvest./2500),
			pagerank(data_net), alpha=0.5,
			ylab=L"\mathrm{share\ PageRank\ (data)}",
			ylabelfontsize=8,
			xticks=(0:8, [L"%$a" for a in 0:8]),
			markerstyle=:star,
			yticks=(0:0.005:0.03, [L"%$a" for a in 0:0.005:0.03]),
			ylim=(0, 0.035),
			xlim=(-0.5, 9)
		),
		layout=(2,1)
	)
	
	postpred_corr = plot(
		left_scatter,
		right_scatter, 
		link=:all,
		bottom_margins=10Plots.mm,
		left_margins=7Plots.mm,
		size=(700, 400), dpi=300
	)

	savefig(postpred_corr, "../images/fig7_post_corr.pdf")

	postpred_corr
end

# ╔═╡ 0c11b8f6-3c9d-4683-a26b-f378d3c7d435
begin
	p01 = scatter(
			log.(1 .+ harvest./2500),
			indegree(data_net),
			ylab=L"\mathrm{indegree}",
			ylabelfontsize=12,
			xticks=false,
			alpha=0.5,
			markercolor=:black,
			yticks=(0:2:14, [L"%$a" for a in 0:2:14]),
			#ylim=(0, 0.035),
			xlim=(-0.5, 9)
		)
	annotate!([1.5], [12], text(L"\mathrm{corr} = %$(round(corspearman(log.(1 .+ harvest./2500), indegree(data_net)), digits=3))", 10))
	
	p02 = scatter(
			log.(1 .+ harvest./2500),
			pagerank(data_net), alpha=0.5,
			ylab=L"\mathrm{share\ PageRank}",
			xlab=L"\log(1 + \mathrm{harvest})",
			ylabelfontsize=12,
			xticks=(0:8, [L"%$a" for a in 0:8]),
			markerstyle=:star,
			yticks=(0:0.005:0.03, [L"%$a" for a in 0:0.005:0.03]),
			ylim=(0, 0.035),
			xlim=(-0.5, 9)
		)
	annotate!([1.5], [0.025], text(L"\mathrm{corr} = %$(round(corspearman(log.(1 .+ harvest./2500), pagerank(data_net)), digits=3))", 10))
	
	corrplot = plot(
		p01,
		p02,
		layout=(2,1)
	)

	#savefig(corrplot, "../images/sup2_corr.pdf")
end

# ╔═╡ 948cf817-93c1-4588-b8d6-77f76e97b85c
md"""
##### Model fit using variance of outdegree
"""

# ╔═╡ b184ffd2-e1c7-495b-b22d-d4f6b6ed742c
begin
	post_pred_var = CSV.read("../data/post_pred_var.csv", DataFrame)
	plot_ppv = plot_post_pred(prior, post_pred_var, data_stats)

	#savefig(plot_ppv, "../images/sup3_post_pred_var.pdf")
end

# ╔═╡ a71a0b62-88cd-4408-8280-742c73d2d8bf
begin
	var1 = "../data/var_fit1.csv"
	var2 = "../data/var_fit2.csv"
	var3 = "../data/var_fit3.csv"
	var4 = "../data/var_fit4.csv"
	var5 = "../data/var_fit5.csv"
	
	varpost = aggregate_posts([var1, var2, var3, var4, var5])
	varpost = varpost[varpost.round .== maximum(varpost.round),:]
	varpost = normalize_pooled_post(varpost)
	CSV.write("../data/var_fit_pooled.csv", varpost)
	
	varpost.c .= varpost.C
	varpost.gamma .= varpost.γ
	varpost.round .= repeat([2], nrow(varpost))
	var_rounds = maximum(varpost.round)
	#varprior = CSV.read("../data/prior_pred.csv", DataFrame)
	
	#each a single‐column plot
	var_pref = plot_posterior_samples(varpost, prior, pref_specs; ncols=1, rounds=rounds, pal=:binary)
	var_comm = plot_posterior_samples(varpost, prior, comm_specs; ncols=1, rounds=rounds, pal=:binary)
	
	#tile them side by side
	post_var = plot(var_pref, var_comm; layout=(1,2), size=(800,600))

	#savefig(post_var, "../images/sup4_post_var.pdf")
end

# ╔═╡ 495fd3a4-d532-42ce-a033-e64edcc3ca62
ess_varpost = ess( varpost.weight ) / nrow(varpost)

# ╔═╡ 4418ffb0-3484-4e20-8e40-97ce3ff3b389
begin
	varpost_end = varpost[varpost.round .== maximum(varpost.round),:]
	
	wv_var = collect(varpost_end.weight); wv ./= sum(wv)
	vμ = [sum(varpost_end[!,c] .* wv_var) for c in param_cols]
	vθs = [
		(
			varpost_end.l[i], varpost_end.C[i], varpost_end.γ[i], varpost_end.beta[i], varpost_end.dens[i], varpost_end.ref[i], varpost_end.strength[i]
		) 
		for i in 1:nrow(varpost_end)
	]
	vi_medoid = argmin([sum((vθs[i] .- vμ).^2) for i in eachindex(vθs)])
	vθ = vθs[vi_medoid]
	
	var_model = initialize_sharing_signals_ywb(
		C = vθ[2], l = vθ[1], γ = vθ[3], beta = vθ[4],
		dens = vθ[5], ref = vθ[6], strength = vθ[7],
		total_ticks = 1000, seed = 7457
	)
	
	run!(
		var_model, 1001, 
		mdata = [
			:l,
			:h,
			:γ,
			:beta,
			:dens,
			:strength,
			:ref,
			:share_density,
			:com_density,
			:share_avg_degree,
			:com_avg_degree,
			:share_avg_clust,
			:com_avg_clust,
			:share_giant_component,
			:com_giant_component,
			:share_diameter,
			:com_diameter,
			:share_avg_pathlength,
			:com_avg_pathlength,
			:com_degree_assort,
			:share_degree_assort,
			:com_wealth_assort,
			:share_wealth_assort,
			:com_endow_deg_corr,
			:share_endow_indeg_corr,
			:com_zero_deg_count,
			:share_zero_outdegree_count,
			:share_zero_indegree_count,
			:com_endow_deg_corr_nonzero,
			:share_endow_deg_corr_nonzero,
			:com_coreness_avg,
			:com_coreness_max,
			:com_coreness_endow_corr,
			:share_coreness_avg,
			:share_coreness_max,
			:share_coreness_endow_corr
	],
		when_model = 1000
	)
	
	nethist_var = net_hist(var_model, data_net, model_inbins=50, model_outbins=50)

	#savefig(nethist_var, "../images/sup5_nethist_var.pdf")
end

# ╔═╡ d7a85cfd-bbeb-40a0-a28f-76e284342306
begin
	netplots_var = plot(
		plot_coreness(var_model.sharenet, (1000 .+ harvest).^1.1, title=L"\mathrm{model}",),
		plot_coreness(data_net, (1000 .+ harvest).^1.1, title=L"\mathrm{data}"),
		layout=(1,2), dpi=300
	)
	
	coreplot_var = plot(
		netplots_var,
		plot_core(data_net, var_model.sharenet, maxcore=11),
		layout=(2,1),
		size=(600,600)
	)

	#savefig(coreplot_var, "../images/sup6_coreplot_var.pdf")
end

# ╔═╡ 65c84d64-3c4c-4e43-8b96-76cf5261c557
begin
	function plot_posterior_pairplot(
		posterior_df::DataFrame,
		prior_df::DataFrame;
		params::Vector{Symbol}=Symbol[],
		xlim_map::Dict{Symbol,Tuple{<:Float64,<:Float64}}=Dict(),
		axis_labels::Dict{Symbol,AbstractString}=Dict(),
		markersize::Number=3,
		alpha::Number=0.4,
		tickfontsize::Int=8,
		xlabelfontsize::Int=10,
		ylabelfontsize::Int=10,
		size::Tuple{Int,Int}=(1500,1500),
		pal::Symbol=:matter,
		dpi::Int=300
	)

	    # 1) find & filter parameters
	    post_cols = Set(Symbol.(names(posterior_df)))
	    if isempty(params)
	        params = sort(collect(post_cols))
	    end
	    good = [p for p in params if p in post_cols]
	    dropped = setdiff(params, good)
	    if !isempty(dropped)
	        @warn "Dropping these (not in posterior): $dropped"
	    end
	    params = good
	    @assert !isempty(params) "No parameters left to plot!"
	
	    n = length(params)
	    grid = Vector{Plots.Plot}(undef, n*n)
		max_round = maximum(posterior_df.round)
	
	    # 2) build each cell
	    for i in 1:n, j in 1:n
	        pi, pj = params[i], params[j]
	        idx    = (i-1)*n + j
	
	        # — base plot —
	        if i == j
	            # diagonal: KDE vs prior
				p = plot(
					plot_param_kde(
						posterior_df,
						prior_df[!, pj],
						pj;
						xlim = get(xlim_map, pj, nothing),
						xlab = "",
						titlefontsize=6,
						rounds = max_round,
						pal = pal
					),
					xticks=false,
					yticks=false
				)
	        else
	            # off‑diagonal: posterior scatter
				post_j = posterior_df[posterior_df.round .== max_round, pj]
				post_i = posterior_df[posterior_df.round .== max_round, pi]
	            p = scatter(
	                post_j,
	                post_i;
					#xlim = (minimum(post_j)-1, maximum(post_j)+1),
					#ylim = (minimum(post_i)-1, maximum(post_i)+1),
	                ms = markersize,
	                alpha = alpha,
	                label = false,
					xticks = false,
					yticks = false,
	            )
	            if haskey(xlim_map, pj)
	                xlims!(p, xlim_map[pj])
	            end
	            if haskey(xlim_map, pi)
	                ylims!(p, xlim_map[pi])
	            end
	        end
	
	        # — axis labels only on edges —
	        if i == n
	            xlabel!(p, get(axis_labels, pj, string(pj)))
	        else
	            xlabel!(p, "")
	            plot!(p, xticks=false)
	        end
	
	        if j == 1
	            ylabel!(p, get(axis_labels, pi, string(pi)))
	        else
	            ylabel!(p, "")
	            plot!(p, yticks=false)
	        end
	
	        # — enforce uniform fonts —
	        plot!(
	            p,
	            xtickfont   = font(tickfontsize),
	            ytickfont   = font(tickfontsize),
	            xguidefont  = font(xlabelfontsize),
	            yguidefont  = font(ylabelfontsize),
	        )
	
	        grid[idx] = p
	    end
	
	    # 3) assemble grid
	    return plot(
	        grid...;
	        layout = (n, n),
	        size   = size,
	        dpi    = dpi,
	        margin = 5Plots.mm,
	        grid   = false,
	    )
	end

	xlims  = Dict(
		:c => (0.0,50.0), 
		:l => (0.0,10.0), 
		:gamma => (0.0,1.0), 
		:beta => (0.0,3.0),
		:dens => (0.0,1.0),
		:ref => (0.0,1.0),
		:dens => (0.0,1.0),
		:strength => (0.0,1.0),
	)
	
	labels = Dict{Symbol, AbstractString}(
		:c => L"\mathrm{c.\ scale\ } (c_\mathrm{eff})",
		:l => L"\mathrm{c.\ elast.\ } (l)",
		:gamma => L"\mathrm{H.\ adv.\ } (γ)",
		:beta  => L"\mathrm{soc.\ reward\ } (\beta)",
		:dens => L"\mathrm{comm.\ dens.\ } (u)",
		:ref => L"\mathrm{central\ H.\ } (\mu)",
		:strength => L"\mathrm{strength\ } (S)"
	)
	
	plot(
		plot_posterior_pairplot(
			post, prior;
			params = [:c, :l, :gamma, :beta, :dens, :ref, :strength],
			xlim_map = xlims,
			axis_labels = labels,
			markersize = 2,
			alpha = 0.5,
			tickfontsize = 8,
			xlabelfontsize = 16,   # ← uniform x-axis label font size
			ylabelfontsize = 16,   # ← uniform y-axis label font size
			size = (1500,1700),
			pal = :binary,
			dpi = 300
		),
		left_margins=10Plots.mm
	)
		
end

# ╔═╡ b802d5cd-534d-41db-9fbe-fb06ff5fba37
# ╠═╡ disabled = true
#=╠═╡
begin
	#each a single‐column plot
	p_pref_unf = plot_posterior_samples(post_unf, prior, pref_specs; ncols=1, rounds=maximum(post_unf.round), pal=:binary)
	p_comm_unf = plot_posterior_samples(post_unf, prior, comm_specs; ncols=1, rounds=maximum(post_unf.round), pal=:binary)
	
	#tile them side by side
	postplot_unf = plot(p_pref_unf, p_comm_unf; layout=(1,2), size=(800,600))
end
  ╠═╡ =#

# ╔═╡ b3837a49-8d56-401d-bf5e-95484812ce49
# ╠═╡ disabled = true
#=╠═╡
begin
	post_pred1.round .= Int.(ones(nrow(post_pred1)))
	prior_pred1.round .= Int.(ones(nrow(prior_pred1)))
	plot_posterior_pairplot(
		prior_pred1, prior_pred1,
		params = [
			:share_avg_degree,
			:share_median_indegree,
	        :share_median_outdegree,
	        :share_avg_clust,
	        :share_coreness_avg,
	        :share_var_indegree,
	        :share_var_outdegree,
	        :share_zero_indegree_count,
	        :pagerank_harvest_corr,
	        :gini_indegree,
		],
		xlim_map = Dict(
			:share_avg_degree => (0.0, 30.0),
			:share_median_indegree => (0.0, 30.0),
	        :share_median_outdegree => (0.0, 30.0),
	        :share_avg_clust => (0.0, 0.3),
	        :share_coreness_avg => (0.0, 15.0),
	        :share_var_indegree => (0.0, 30.0),
	        :share_var_outdegree => (0.0, 30.0),
	        :share_zero_indegree_count => (0.0, 30.0),
	        :pagerank_harvest_corr => (-1.0, 1.0),
	        :gini_indegree => (0.0, 1.0),
		),
		axis_labels = labels,
		markersize = 2,
		alpha = 0.5,
		tickfontsize = 8,
		xlabelfontsize = 16,   # ← uniform x-axis label font size
		ylabelfontsize = 16,   # ← uniform y-axis label font size
		size = (1500,1700),
		pal = :binary,
		dpi = 300
	)
end
  ╠═╡ =#

# ╔═╡ fc3e8e24-4276-4724-b4a0-5f1b68b543f5
# ╠═╡ disabled = true
#=╠═╡
begin
	function ABC(sim_stats, data_stats, col_list)
	    # Remove rows with any NaN values in the selected columns
	    sim_stats = filter(row -> all(!isnan, row[col_list]), sim_stats)
	    
	    # Compute the covariance matrix of the simulated summary statistics and its inverse
	    cov_matrix = cov(Matrix(sim_stats[:, col_list]))
	    inv_cov_matrix = inv(cov_matrix)
	    
	    # Compute the Mahalanobis distance for each simulated row relative to every observed row,
	    # then take the mean distance as the distance metric.
	    sim_stats.mah_dist = [
	        mean([
	            sqrt(
	                (Matrix(sim_stats[i:i, col_list])[:] - Matrix(data_stats[j:j, col_list])[:])' *
	                inv_cov_matrix *
	                (Matrix(sim_stats[i:i, col_list])[:] - Matrix(data_stats[j:j, col_list])[:])
	            ) for j in 1:nrow(data_stats)
	        ])
	        for i in 1:nrow(sim_stats)
	    ]
	    
	    # Set the threshold as the 0.1th percentile of the Mahalanobis distances.
	    ε_mah = quantile(sim_stats.mah_dist, 0.0005)
	    
	    # Filter simulated results that are within the threshold distance
	    sim_stats_filtered = sim_stats[sim_stats.mah_dist .< ε_mah, :]
		sim_stats_filtered.round = repeat([2], nrow(sim_stats_filtered))
	    
	    return sim_stats_filtered
	end

	posterior = ABC(prior1, data_stats, col_list)
	posterior.γ = posterior.gamma
	
	# make each a single‐column plot:
	p_pref = plot_posterior_samples(posterior, prior1, pref_specs; ncols=1, rounds=2)
	p_comm = plot_posterior_samples(posterior, prior1, comm_specs; ncols=1, rounds=2)
	
	# then tile them side by side:
	plot(p_pref, p_comm; layout=(1,2), size=(800,600))
end
  ╠═╡ =#

# ╔═╡ 1f4b3687-8d40-4eda-817e-52b6d287d406
# ╠═╡ disabled = true
#=╠═╡
plot(
	plot_posterior_pairplot(
	  posterior, prior;
	  params          = [:C, :l, :γ, :beta, :dens, :ref, :strength],
	  xlim_map        = xlims,
	  axis_labels     = labels,
	  markersize      = 2,
	  alpha           = 0.5,
	  tickfontsize    = 8,
	  xlabelfontsize  = 16,   # ← uniform x-axis label font size
	  ylabelfontsize  = 16,   # ← uniform y-axis label font size
	  size            = (1500,1700),
	  dpi             = 300
	),
	left_margins=10Plots.mm
)
  ╠═╡ =#

# ╔═╡ 3a0dacc8-b717-463c-9a38-1923e4d3803c
# ╠═╡ disabled = true
#=╠═╡
plot(
	plot_prior_post(
		prior_pred0, 
		post_pred0, 
		data_stats,
		col_list1, 
		ylim = (-0.1, 90),
		xnames = ([1,2,3,4,5],[L"\mathrm{MEAN.\ DEG}", L"\mathrm{MED.\ IN}", L"\mathrm{MED.\ OUT}", L"\mathrm{CORE}", L"\mathrm{ZERO}"]),
		ynames = ([0, 10, 20, 30, 40, 50, 60, 70, 80], [L"%$a" for a in [0, 10, 20, 30, 40, 50, 60, 70, 80]]),
	),
	plot_prior_post(
		prior_pred0, 
		post_pred0,
		data_stats,
		col_list2, 
		xnames = ([1], [L"\mathrm{CLUST}"]),
		ynames = ([0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15], [L"%$a" for a in [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15]]),
		ylabel = "", 
		offset = 0.35, 
		xlim = (-0.1, 2), 
		ylim = (-0.0003, 0.17),
		legend = false
	),
	plot_prior_post(
		prior_pred0, 
		post_pred0,
		data_stats,
		col_list3, 
		xnames = (
		[1],[L"\mathrm{VAR.\ OUT}"]
		),
		ynames = ([0, 25, 50, 75, 100, 125, 150, 175], [L"%$a" for a in [0, 25, 50, 75, 100, 125, 150, 175]]),
		ylabel="", 
		offset=0.2, 
		xlim=(0.5, 1.5), 
		ylim = (-0.5, 200),
		legend=false
	),
	layout = @layout([a{0.6w} b{0.2w} c{0.2w}]), bottom_margins=8Plots.mm, left_margins=2Plots.mm,size=(800, 400)
)
  ╠═╡ =#

# ╔═╡ f7412c10-12fb-4e20-973f-9322253e9975
# ╠═╡ disabled = true
#=╠═╡
begin
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
	
	function ref_summaries(sharenet::DiGraph, harvest::Vector{<:Real})
		gR = sharenet
		auth  = pagerank(gR)                       
		core  = core_number(gR)
		k_in  = indegree(sharenet)
	
		ρ_auth     = corspearman(auth, harvest)
		ρ_core_in  = corspearman(core, harvest)
		gini_in    = gini(float.(k_in))
	
		return (; ρ_auth, ρ_core_in, gini_in)      
	end
	
	ref_summaries(data_net, harvest)
end
  ╠═╡ =#

# ╔═╡ a54a6c23-777c-49d7-be05-78584ea936a2
# ╠═╡ disabled = true
#=╠═╡
begin
	function find_nan_locations(df::DataFrame)
	    nan_locations = []
	    for col in names(df)
	        for i in 1:size(df, 1)
	            if isnan(df[i, col])
	                push!(nan_locations, (i, col))
	            end
	        end
	    end
	    return nan_locations
	end
	
	find_nan_locations(prior1)
end
  ╠═╡ =#

# ╔═╡ b54ffd56-09eb-47ee-bce6-10fdd265e17c
# ╠═╡ disabled = true
#=╠═╡
begin
	data_net, harvest = datanet()

	prior1 = CSV.read("../data/analysis_1.csv", DataFrame)
	posterior1 = CSV.read("../data/analysis_1_filtered.csv", DataFrame)

	prior1 = prior1[prior1.h .< 20 .&& prior1.l .< 20 .&& prior1.gamma .< 20 .&& prior1.beta .< 20, :]
	
	function binned_summary_plot(x::AbstractVector, y::AbstractVector;
	                             m::Int=50,
	                             quantiles::Tuple{<:Real,<:Real}=(0.25,0.75),
	                             xlabel::AbstractString="",
	                             ylabel::AbstractString="",
	                             title::String="",
	                             kwargs...)
	
	    @assert length(x) == length(y) "x and y must have same length"
	
	    # 1) bin edges & centers
	    xmin, xmax = minimum(x), maximum(x)
	    edges = range(xmin, stop=xmax, length=m+1)
	    centers = [(edges[i] + edges[i+1]) / 2 for i in 1:m]
	
	    medians = Float64[]
	    lowers  = Float64[]
	    uppers  = Float64[]
	
	    # 2) compute statistics in each bin
	    qlow, qhigh = quantiles
	    for i in 1:m
	        lo, hi = edges[i], edges[i+1]
	        inds = findall(t -> (t ≥ lo) && (t < hi), x)
	        if !isempty(inds)
	            ys = y[inds]
	            push!(medians, median(ys))
	            qs = quantile(ys, [qlow, qhigh])
	            push!(lowers, qs[1])
	            push!(uppers, qs[2])
	        else
	            # if a bin is empty, push NaN so plotting skips it
	            push!(medians, NaN)
	            push!(lowers,  NaN)
	            push!(uppers,  NaN)
	        end
	    end
	
	    # 3) make the plot: median line + ribbon
	    plot(
	      centers, medians;
	      ribbon    = (medians .- lowers, uppers .- medians),
		  legend 	= false,
	      label     = "median ± quantiles",
	      xlabel    = xlabel,
	      ylabel    = ylabel,
	      title     = title,
		  color 	= :black,
		  lw 		= 2,
		  grid 		= false,
	      kwargs...
	    )
	end
	
	plot(
		plot(	
			binned_summary_plot(
				prior1.h, prior1.share_avg_degree, 
				ylabel=L"\mathrm{mean\ degree}", 
				xlim=(0,20), 
				xticks=false, yticks=(0:20:60, [L"%$a" for a in 0:20:60])
			),
			binned_summary_plot(
				prior1.l, prior1.share_avg_degree, 
				ylabel="", 
				xlim=(0,20), yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.gamma, prior1.share_avg_degree, 
				ylabel="", 
				xlim=(0,20), yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.beta, prior1.share_avg_degree, 
				ylabel="", 
				xlim=(0,20), yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.dens, prior1.share_avg_degree, 
				ylabel="", yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.ref, 
				prior1.share_avg_degree, 
				ylabel="", yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.strength, prior1.share_avg_degree, 
				ylabel="", yticks=false, xticks=false
			),
			layout=(1,7), link=:y,
		),
		plot(	
			binned_summary_plot(
				prior1.h, prior1.share_avg_clust, 
				ylabel=L"\mathrm{clustering}", xlim=(0,20), 
				xticks=false, yticks=(0.0:0.05:0.2, [L"%$a" for a in 0.0:0.05:0.2])
			),
			binned_summary_plot(
				prior1.l, prior1.share_avg_clust, 
				ylabel="", xlim=(0,20), 
				yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.gamma, prior1.share_avg_clust, 
				ylabel="", xlim=(0,20), 
				yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.beta, prior1.share_avg_clust, 
				ylabel="", xlim=(0,20), 
				yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.dens, prior1.share_avg_clust, 
				ylabel="", yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.ref, prior1.share_avg_clust, 
				ylabel="", yticks=false, xticks=false
			),
			binned_summary_plot(
				prior1.strength, prior1.share_avg_clust, 
				ylabel="", yticks=false, xticks=false
			),
			layout=(1,7), link=:y,
		),
		plot(	
			binned_summary_plot(
				prior1.h, prior1.share_coreness_avg, 
				ylabel=L"\mathrm{coreness}", xlabel=L"\mathrm{H. elasticity\ } (h)",
				xlabelfontsize=8,
				xlim=(0,20), yticks=(0:25:75, [L"%$a" for a in 0:25:75]),
				xticks=(0:5:20, [L"%$a" for a in 0:5:20])
			),
			binned_summary_plot(
				prior1.l, prior1.share_coreness_avg, 
				ylabel="", xlabel=L"\mathrm{cost elasticity\ } (l)",
				xlabelfontsize=8,
				xlim=(0,20), yticks=false,
				xticks=(0:5:20, [L"%$a" for a in 0:5:20])
			),
			binned_summary_plot(
				prior1.gamma, prior1.share_coreness_avg, 
				ylabel="", xlabel=L"\mathrm{H. advantage\ } (\gamma)",
				xlabelfontsize=8,
				xlim=(0,20), yticks=false,
				xticks=(0:5:20, [L"%$a" for a in 0:5:20])
			),
			binned_summary_plot(
				prior1.beta, prior1.share_coreness_avg, 
				ylabel="", xlabel=L"\mathrm{social\ reward\ } (\beta)",
				xlabelfontsize=8,
				xlim=(0,20), yticks=false,
				xticks=(0:5:20, [L"%$a" for a in 0:5:20])
			),
			binned_summary_plot(
				prior1.dens, prior1.share_coreness_avg, 
				ylabel="", yticks=false, xlabel=L"\mathrm{density\ } (u)",
				xlabelfontsize=8,
				xticks=(0.0:0.5:1.0, [L"%$a" for a in 0.0:0.5:1.0])
			),
			binned_summary_plot(
				prior1.ref, prior1.share_coreness_avg, 
				ylabel="", yticks=false, xlabel=L"\mathrm{central\ harvest\ } (\mu)",
				xlabelfontsize=8,
				xticks=(0.0:0.5:1.0, [L"%$a" for a in 0.0:0.5:1.0])
			),
			binned_summary_plot(
				prior1.strength, prior1.share_coreness_avg, 
				ylabel="", yticks=false, xlabel=L"\mathrm{strength\ } (S)",
				xlabelfontsize=8,
				xticks=(0.0:0.5:1.0, [L"%$a" for a in 0.0:0.5:1.0])
			),
			layout=(1,7), link=:y,
		),
		layout=(3,1), size=(800, 600)
	)
end
  ╠═╡ =#

# ╔═╡ f05af632-7a78-49fd-8b4f-a90c1ebf3d2b
# ╠═╡ disabled = true
#=╠═╡
begin
	trimmed_prior = prior1[prior1.dens .< 0.15, :]
	
	plot(
		binned_summary_plot(
			trimmed_prior.dens, trimmed_prior.share_zero_indegree_count, 
			ylabel=L"\mathrm{zero\ \ indegree\ \ count}", yticks=(0:25:100, [L"%$a" for a in 0:25:100]),
			xticks=(0.0:0.05:0.15, [L"%$a" for a in 0.0:0.05:0.15]),
			xlabel=L"\mathrm{density\ } (u)"
		),
		binned_summary_plot(
			trimmed_prior.ref, trimmed_prior.share_zero_indegree_count, 
			#ylabel=L"\mathrm{zero\ \ indegree\ \ count}",
			yticks=(0:10:50, [L"%$a" for a in 0:10:50]),
			xticks=(0.0:0.25:1.0, [L"%$a" for a in 0.0:0.25:1.0]),
			xlabel=L"\mathrm{central\ harvest\ } (\mu)"
		),
		binned_summary_plot(
			trimmed_prior.strength, trimmed_prior.share_zero_indegree_count, 
			ylabel="", yticks=false,
			xticks=(0.0:0.25:1.0, [L"%$a" for a in 0.0:0.25:1.0]),
			xlabel=L"\mathrm{strength\ } (S)"
		),
		layout=(1,3), link=:y, size=(600, 300), plot_title=L"\mathrm{zero\ indegree\ count\ at\ } u < 0.15"
	)
end
  ╠═╡ =#

# ╔═╡ 5d9244fc-cb9e-4214-b50f-0df5535186f6
# ╠═╡ disabled = true
#=╠═╡
begin
	function run_model_expt(
		params::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,Int})
	    l_val, h_val, gamma_val, beta_val, dens_val, ref_val, strength_val, seed_val = params
	        tmax = 1000  # total simulation ticks (adjust as needed)
	        
	        # Initialize the model (replace with your actual function):
	        model = initialize_sharing_signals_lite(
	            l = l_val,
	            h = h_val,
	            γ = gamma_val,
	            beta = beta_val,
	            dens = dens_val,
	            ref = ref_val,
	            strength = strength_val,
	            seed = seed_val,
	            total_ticks = tmax
	        )
	        
	        # Run the simulation
	        for _ in 1:tmax
	            step!(model)
	        end
	        
	        # Collect model outputs
	        p = model
	        return (
	            l = l_val,
	            h = h_val,
	            gamma = gamma_val,
	            beta = beta_val,
	            dens = dens_val,
	            ref = ref_val,
	            strength = strength_val,
	            seed = seed_val,
	            final_tick = p.tick,
	            com_density = p.com_density,
	            com_avg_degree = p.com_avg_degree,
	            com_avg_clust = p.com_avg_clust,
	            com_giant_component = p.com_giant_component,
	            com_coreness_avg = p.com_coreness_avg,
	            com_coreness_max = p.com_coreness_max,
	            com_coreness_endow_corr = p.com_coreness_endow_corr,
	            com_endow_deg_corr_nonzero = p.com_endow_deg_corr_nonzero,
	            com_endow_deg_corr = p.com_endow_deg_corr,
	            com_wealth_assort = p.com_wealth_assort,
	            com_degree_assort = p.com_degree_assort,
	            com_zero_deg_count = p.com_zero_deg_count,
	            com_avg_pathlength = p.com_avg_pathlength,
	            com_diameter = p.com_diameter,
	            share_density = p.share_density,
	            share_avg_degree = p.share_avg_degree,
	            share_median_indegree = p.share_median_indegree,
	            share_median_outdegree = p.share_median_outdegree,
	            share_var_degree = p.share_var_degree,
	            share_var_indegree = p.share_var_indegree,
	            share_var_outdegree = p.share_var_outdegree,
	            share_avg_clust = p.share_avg_clust,
	            share_giant_component = p.share_giant_component,
	            share_coreness_avg = p.share_coreness_avg,
	            share_coreness_max = p.share_coreness_max,
	            share_coreness_endow_corr = p.share_coreness_endow_corr,
	            share_endow_deg_corr_nonzero = p.share_endow_deg_corr_nonzero,
	            share_endow_deg_corr = p.share_endow_deg_corr,
	            share_degree_assort = p.share_degree_assort,
	            share_wealth_assort = p.share_wealth_assort,
	            share_zero_outdegree_count = p.share_zero_outdegree_count,
	            share_zero_indegree_count = p.share_zero_indegree_count,
	            share_avg_pathlength = p.share_avg_pathlength,
	            share_diameter = p.share_diameter,
	            total_warns = p.total_warns
	        )
	end
	n_param = 10
	parameters = ( #ALTER THIS DICTIONARY TO DEFINE PARAMETER DISTRIBUTIONS
	    repeat([3.0], n_param),
		repeat([3.0], n_param),
		repeat([0.5], n_param),
		repeat([0.5], n_param),
		repeat([0.1], n_param),
		repeat([0.35], n_param),
		repeat([0.95], n_param),
		rand(1:100000, n_param)
	)
	parameters = collect(zip(parameters...))
	test_data = DataFrame(pmap(run_model_expt, parameters))

	test_fit = ABC(sim_stats, test_data, col_list)
	plot_posterior_samples(test_fit)
end
  ╠═╡ =#

# ╔═╡ d3f36c06-a8d9-4a0e-a39c-c1f94bd1e171
# ╠═╡ disabled = true
#=╠═╡
plot(
	0:10,
	1 .- exp.(-1.18 .* (0:10|>collect)),
	legend=false
)
  ╠═╡ =#

# ╔═╡ b0ce249a-c5c2-4276-96b1-3d748ab29382
# ╠═╡ disabled = true
#=╠═╡
begin
	modelplot = plot(
		scatter(
			log.(1 .+ (harvest./2000)),
			[a.indeg for a in allagents(model)|>collect],
			legend=false,
			xlab=L"\mathrm{log\ endowment}",
			ylab=L"\mathrm{indegree}",
			alpha=0.5,
			ylim=(0, 17),
			color=palette(:tokyo10)[3],
			grid=false,
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
			yticks=(0:5:15, [L"%$a" for a in 0:5:15]),
		),
		plot(
			histogram(
				[a.indeg for a in allagents(model)|>collect], xlab=L"\mathrm{indegree}", 
				legend=false, 
				xlim=(0, 17),
				color=palette(:tokyo10)[3], 
				bins=65,
				grid=false,
				xticks=(0:5:15, [L"%$a" for a in 0:5:15]),
				yticks=(0:5:20, [L"%$a" for a in 0:5:20]),
			),
			histogram(
				[a.outdeg for a in allagents(model)|>collect], xlab=L"\mathrm{outdegree}", 
				legend=false, xlim=(0, 30), 
				ylim=(0, 35), bins=25, 
				color=palette(:tokyo10)[3],
				grid=false,
				xticks=(0:10:30, [L"%$a" for a in 0:10:30]),
				yticks=(0:10:30, [L"%$a" for a in 0:5:35]),
			),
			layout=(2,1)
		)
	)
	
	dataplot = plot(
		scatter(
			log.(1 .+ (harvest./2000)),
			indegree(data_net),
			bins=30, 
			xlab=L"\mathrm{log\ endowment}",
			ylab=L"\mathrm{indegree}",
			alpha=0.5,
			legend=false,
			color=palette(:tokyo10)[7],
			grid=false,
			xticks=(0:2:8, [L"%$a" for a in 0:2:8]),
			yticks=(0:5:15, [L"%$a" for a in 0:5:15]),
		),
		plot(
			histogram(
				indegree(data_net), 
				xlab=L"\mathrm{indegree}", 
				legend=false, bins=20, 
				color=palette(:tokyo10)[7],
				grid=false,
				xticks=(0:5:15, [L"%$a" for a in 0:5:15]),
				yticks=(0:5:20, [L"%$a" for a in 0:5:20]),
			),
			histogram(
				outdegree(data_net), 
				xlab=L"\mathrm{outdegree}", 
				legend=false, bins=30, 
				ylim=(0,35), 
				color=palette(:tokyo10)[7],
				grid=false,
				xticks=(0:10:30, [L"%$a" for a in 0:10:30]),
				yticks=(0:10:30, [L"%$a" for a in 0:5:35]),
			),
			layout=(2,1)
		)
	)
	
	plot(
		plot(modelplot, plot_title=L"\mathrm{model\ network}", margins=2Plots.mm),
		plot(dataplot, plot_title=L"\mathrm{data\ network}", margins=2Plots.mm),
		size=(700, 400)
	)
end
  ╠═╡ =#

# ╔═╡ a7ab5745-fabc-4164-b720-23d178ca3456
# ╠═╡ disabled = true
#=╠═╡
begin
	coreness = []
		for i in 0:15
			push!(coreness, length(k_core(model.sharenet, i)) - length(k_core(model.sharenet, i+1)))
		end
	plot(0:15, coreness, lw=2, legend=false, xlabel="k-shell")

	corona = []
	for i in 0:15
		push!(corona, length(k_shell(model.sharenet, i)))
	end
	#sim_coreplot = plot(1:20, corona, lw=2, legend=false, xlabel="k-core", color="dark red")
	plot!(0:15, corona, lw=2, legend=false, xlabel="k-shell", color="dark red")
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─887a17bf-76fd-445b-a82d-ff5c5e45d3f9
# ╟─a5a9ecc8-5750-40fe-b433-e7f21d5b8fc8
# ╟─da8d2640-88d5-406a-bf74-32575c229963
# ╟─13a5a8bc-e920-46c8-b24e-0fee93f4cf8b
# ╟─d4bb1701-efa7-4bb0-ab69-73e0b7b8f5ac
# ╟─e22c7a83-bff9-4766-bcf4-071077150b13
# ╟─56af3b89-c9aa-405e-9deb-cabb77e69993
# ╟─cc9089be-7dc3-400a-91dc-1f5c442ee41d
# ╟─988cd7d4-b6f4-4579-bdca-45ee88854f38
# ╟─b2d0bb8c-00a8-4930-a7d2-4ccf24ca752d
# ╟─88f07c6d-ad7f-4d9f-845b-3ef91cc59f28
# ╟─a657563c-fa1e-455d-ac70-3fabdd52281c
# ╟─88fefd48-f055-4a0e-b56f-1f38d8eae9bc
# ╟─2c6e7505-8036-4a84-bad8-19d861034a53
# ╟─28f2262c-56e8-429c-81e9-ecdc432cd691
# ╟─d8a0ccfc-c28f-446c-b61d-e0ca8933dec8
# ╟─85943361-ab1e-4266-9c59-c4184b49948f
# ╟─8db4d3f1-5060-46a6-a9ff-4ca2941b5247
# ╟─0c11b8f6-3c9d-4683-a26b-f378d3c7d435
# ╟─948cf817-93c1-4588-b8d6-77f76e97b85c
# ╟─b184ffd2-e1c7-495b-b22d-d4f6b6ed742c
# ╟─a71a0b62-88cd-4408-8280-742c73d2d8bf
# ╟─495fd3a4-d532-42ce-a033-e64edcc3ca62
# ╟─4418ffb0-3484-4e20-8e40-97ce3ff3b389
# ╟─d7a85cfd-bbeb-40a0-a28f-76e284342306
# ╟─65c84d64-3c4c-4e43-8b96-76cf5261c557
# ╟─b802d5cd-534d-41db-9fbe-fb06ff5fba37
# ╟─b3837a49-8d56-401d-bf5e-95484812ce49
# ╟─fc3e8e24-4276-4724-b4a0-5f1b68b543f5
# ╟─1f4b3687-8d40-4eda-817e-52b6d287d406
# ╟─3a0dacc8-b717-463c-9a38-1923e4d3803c
# ╟─f7412c10-12fb-4e20-973f-9322253e9975
# ╟─a54a6c23-777c-49d7-be05-78584ea936a2
# ╟─b54ffd56-09eb-47ee-bce6-10fdd265e17c
# ╟─f05af632-7a78-49fd-8b4f-a90c1ebf3d2b
# ╟─5d9244fc-cb9e-4214-b50f-0df5535186f6
# ╟─d3f36c06-a8d9-4a0e-a39c-c1f94bd1e171
# ╟─b0ce249a-c5c2-4276-96b1-3d748ab29382
# ╟─a7ab5745-fabc-4164-b720-23d178ca3456
