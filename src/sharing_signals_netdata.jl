using Graphs, DataFrames, CSV, LinearAlgebra, SparseArrays

function construct_network(n, edgelist; delete_vertices = false, ids_delete=[0], node_color=1)
    eT = []
    GAB = edgelist.GAB
    GBA = edgelist.GBA
    A = edgelist.A
    B = edgelist.B
    for i in 1:nrow(edgelist)
        if GAB[i] == 1
            push!(eT, (A[i], B[i]))
        end
        if GBA[i] == 1
            push!(eT, (B[i], A[i]))
        end
    end

    eT = Edge.(eT)
    
    gT = SimpleDiGraph(n)

    for e in eT
        add_edge!(gT, e)
    end

    if delete_vertices
        rem_vertices!(gT, ids_delete, keep_order=true)
    end

    return [gT, eT]
end

function datanet()
	network_data = CSV.read("../netdata/cf_sharing_edgelist.csv", DataFrame)
	ses_attributes = CSV.read("../netdata/ses_variables.csv", DataFrame)

	network_data_filtered = filter(a -> a.A_to_B_rA != "NA", network_data)
	network_data_filtered = filter(a -> a.A_to_B_rB != "NA", network_data_filtered)
	network_data_filtered = filter(a -> a.B_to_A_rA != "NA", network_data_filtered)
	network_data_filtered = filter(a -> a.B_to_A_rB != "NA", network_data_filtered)
	
	for c ∈ [
		"A_to_B_rA",
		"A_to_B_rB",
		"B_to_A_rA",
		"B_to_A_rB"
	]
    	network_data_filtered[!, c]= parse.(Int64, network_data_filtered[!, c])
	end

	network_data_filtered[!, "GAB"] = network_data_filtered[!, "A_to_B_rA"] .+ network_data_filtered[!, "A_to_B_rB"]
	network_data_filtered[!, "GBA"] = network_data_filtered[!, "B_to_A_rA"] .+ network_data_filtered[!, "B_to_A_rB"]

	discount(n) = n == 2 ? 1 : n
	
	network_data_filtered[!,"GAB"] = discount.(network_data_filtered.GAB)
	network_data_filtered[!,"GBA"] = discount.(network_data_filtered.GBA)

	network_data_f = deepcopy(network_data_filtered)

	extracted_ids = filter(i -> i.harvest != "NA", ses_attributes)
	extracted_ids = extracted_ids.HID
	net_filter = setdiff(1:179, extracted_ids)

	#zero_harvest = filter(i -> i.harvest == 0, ses_attributes)
	#zero_ids = zero_harvest.HID
	
	#maskA2 = in.(network_data_f.A, Ref(zero_harvest[:,:HID]))
	#maskB2 = in.(network_data_f.B, Ref(zero_harvest[:,:HID]))
	
	#network_data_f[maskA2, :GAB] .= 0
	#network_data_f[maskB2, :GBA] .= 0

	surinv = CSV.read("../netdata/Survey_individuals.csv", DataFrame)
	for col in names(surinv)
		surinv[:,col] = replace(surinv[:,col], "NULL" => "0")
	end
	for col in names(surinv)[26:30]
	    surinv[!, col] = string.(surinv[!, col])
	    surinv[!, col] = parse.(Int, surinv[!, col])
	end
	surinv_comb = combine(
		groupby(surinv, [:HID_fk]),
		:HID_fk => mean => :HID,
		#:Fishes => sum => :fish, #we don't include fish
		:Seal_hunting => sum => :seal,
		:Beluga_hunting => sum => :beluga,
		:Caribou_hunting => sum => :caribou,
		:Bird_hunting => sum => :bird,
	)
	surinv_comb.total = [ sum(row) for row in eachrow(surinv_comb[:, 3:end]) ]
	surinv_comb[surinv_comb.total .== 0, :HID_fk]

	maskA = in.(network_data_f.A, Ref(surinv_comb[surinv_comb.total .== 0, :HID_fk]))
	maskB = in.(network_data_f.B, Ref(surinv_comb[surinv_comb.total .== 0, :HID_fk]))

	network_data_f[maskA, :GAB] .= 0
	network_data_f[maskB, :GBA] .= 0


	data_net = construct_network(
		179,
		network_data_f,
		delete_vertices=true,
		ids_delete=net_filter
	)

	data_net_unfiltered = construct_network(
		179,
		network_data_filtered,
		delete_vertices=true,
		ids_delete=net_filter
	)

	harvest = ses_attributes.harvest

	return data_net[1], harvest, data_net_unfiltered[1]
end

function extract_netdata(g, harvest)
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

	colnames = [
		:share_avg_degree, 
		:share_avg_clust, 
		:share_median_outdegree, 
		:share_var_outdegree, 
		:share_maximum_indegree, 
		:share_zero_indegree_count, 
		:pagerank_harvest_corr, 
		:gini_indegree
	]

	data = [
		[mean(degree(g))],
		[global_clustering_coefficient(g)],
		[median(outdegree(g))],
		[var(outdegree(g))],
		[maximum(indegree(g))],
		[length(filter(x -> x == 0, indegree(g)))],
		[corspearman(pagerank(g), harvest)],
		[gini(indegree(g))]
	]

	return DataFrame(data, colnames)
end

function reciprocity_igraph2(g; mode::Symbol = :dyad)
	# sparse adjacency matrix of g (out-edges => rows are sources, cols are targets)
	A  = Graphs.adjacency_matrix(g; dir = :out)
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

function gini2(v::AbstractVector{<:Real})
	n = length(v)
	n == 0 && return 0.0
	s = sort(float.(v))                # ascending
	total = sum(s)
	total == 0 && return 0.0
	coef = 2 / (n * total)
	g = coef * sum( j * s[j] for j in 1:n ) - (n + 1) / n
	return g                           # 0  ≤ g ≤ 1
end

function extract_netdata_full(g, harvest)
	

	colnames = [
		:share_avg_degree, 
		:share_avg_clust, 
		:share_median_indegree,
		:share_median_outdegree, 
		:share_var_indegree,
		:share_mad_indegree,
		:share_iqr_indegree,
		:share_var_outdegree,
		:share_mad_outdegree,
		:share_iqr_outdegree, 
		:share_maximum_indegree, 
		:share_zero_indegree_count,
		:share_coreness_avg, 
		:pagerank_harvest_corr, 
		:gini_indegree,
		:reciprocity_dyad,
		:reciprocity_edge
	]

	data = [
		[mean(degree(g))],
		[global_clustering_coefficient(g)],
		[median(indegree(g))],
		[median(outdegree(g))],
		[var(indegree(g))],
		[mad(indegree(g))],
		[iqr(indegree(g))],
		[var(outdegree(g))],
		[mad(outdegree(g))],
		[iqr(outdegree(g))],
		[maximum(indegree(g))],
		[length(filter(x -> x == 0, indegree(g)))],
		[mean(core_number(g))],
		[corspearman(pagerank(g), harvest)],
		[gini2(indegree(g))],
		[reciprocity_igraph2(g)],
		[reciprocity_igraph2(g, mode=:edge)],
	]

	return DataFrame(data, colnames)
end