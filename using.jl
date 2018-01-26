# using Lint # check code
# using JLD # save variables
# using OptimizationProblems
using NLPModels
using LinearOperators
using Plots # graph
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
loggerCUTEst = get_logger("launch_CUTEst_pbs")
loggerCUTEstLS = get_logger("launch_CUTEstLS_pbs")
loggerind = get_logger("launch_ind_pbs")
loggercvx = get_logger("launch_cvx_pbs")
loggerNLS = get_logger("NLSProblems")

loggerCR.level = MiniLogging.ERROR
loggerTR.level = MiniLogging.ERROR
loggerCG.level = MiniLogging.ERROR
# loggerminres.level = MiniLogging.ERROR
loggerTRLS.level = MiniLogging.ERROR

include("CG.jl")
include("CR.jl")
# include("minres.jl")
include("TrustRegion.jl")
include("TrustRegionLS.jl")
include("TRCR.jl")
include("TRCG.jl")
# include("TRminres.jl")
include("TRCGLS.jl")
include("TRCRLS.jl")
include("TRLSQR.jl")
include("TRLSMR.jl")
