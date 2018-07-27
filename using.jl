# using Lint # check code
# using JLD # save variables
using NLPModels
using LinearOperators
using Plots # graph
Plots.pyplot()
using BenchmarkTools
using BenchmarkProfiles
using MiniLogging

# basic configuration. The root logger level is then INFO
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")
loggerCR = get_logger("CR")
loggerCG = get_logger("CG")
# loggerminres = get_logger("minres")
loggerTR = get_logger("TrustRegion")
loggerTRLS = get_logger("TrustRegionLS")
loggerTRCUTEst = get_logger("TrustRegionCUTEst")
# loggerTRts = get_logger("TrustRegion_two_solvers")

# loggerCUTEst = get_logger("launch_CUTEst_pbs")
# loggerCUTEstCG = get_logger("launch_CUTEst_CG")
# loggerCUTEstCR = get_logger("launch_CUTEst_CR")
# loggerCUTEstLS = get_logger("launch_CUTEstLS_pbs")
loggerind = get_logger("launch_ind_pbs")
loggercvx = get_logger("launch_cvx_pbs")
# loggerNLS = get_logger("NLSProblems")


loggerTR.level = MiniLogging.ERROR
# loggerTRCUTEst.level = MiniLogging.ERROR
# loggerTRts.level = MiniLogging.ERROR
loggerCG.level = MiniLogging.ERROR
# loggerminres.level = MiniLogging.ERROR
loggerTRLS.level = MiniLogging.ERROR

include("CG.jl")
include("CR.jl")
# include("minres.jl")
include("TrustRegion.jl")
include("TrustRegionCUTEst.jl")
# include("TrustRegion_two_solvers.jl")
include("TrustRegionLS.jl")
include("TRCR.jl")
include("TRCG.jl")
# include("TRCG_CUTEst.jl")
# include("TRCR_CUTEst.jl")
# include("TRminres.jl")
# include("TRCGLS.jl")
# include("TRCRLS.jl")
include("TRLSQR.jl")
include("TRLSMR.jl")
