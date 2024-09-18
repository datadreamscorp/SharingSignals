#############
#ANALYSIS 0
@everywhere using Pkg
@everywhere Pkg.activate("..")
@everywhere using SharingSignals, Agents, Statistics, CSV

@everywhere total_ticks = 1000

@everywhere begin #INCLUDE MODEL CODE AND NECESSARY LIBRARIES

	parameters = Dict( #ALTER THIS DICTIONARY TO DEFINE PARAMETER DISTRIBUTIONS
    :N => [50, 100],
    :B => 1.0,
    :C => 0.0:0.1:1.0|>collect,
    :sigma => 1.0:0.1:5.0|>collect,
    :dens => 0.0:0.05:1.0|>collect,
    :beta => 0.0:0.5:10.0|>collect,
    :seed => 1000:1010|>collect
)

	adata = [
        :B,
        :payoff,
        :total_payoff,
        :comdeg,
        :outdeg,
        :indeg
        ]

end

#USE THIS LINE AFTER DEFINITIONS TO BEGIN PARAMETER SCANNING
_, mdf = paramscan(
            parameters, initialize_sharing_signals;
            adata=adata,
            n = total_ticks,
			parallel=true,
			when_model = total_ticks,
			showprogress = true
	)

CSV.write("../data/analysis_0.csv", mdf)
