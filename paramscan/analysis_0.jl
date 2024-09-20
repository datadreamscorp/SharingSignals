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
    :C => 0.1:0.1:1.0|>collect,
    :sigma => 1.0:0.25:3.0|>collect,
    :dens => 0.0:0.05:0.5|>collect,
    :beta => 0.5:0.5:3.0|>collect,
    :seed => 1000:1010|>collect
)

	adata = [
        :endow,
        :payoff,
        :total_payoff,
        :comdeg,
        :outdeg,
        :indeg
        ]

end

#USE THIS LINE AFTER DEFINITIONS TO BEGIN PARAMETER SCANNING
adf, _ = paramscan(
            	parameters, initialize_sharing_signals;
            	adata = adata,
            	n = total_ticks,
		parallel=true,
		when = [total_ticks],
		showprogress = true
	)

CSV.write("../data/analysis_0.csv", adf)
