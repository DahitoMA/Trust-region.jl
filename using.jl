# using Lint # check code
# using JLD # save variables
using OptimizationProblems
using NLPModels
using LinearOperators
# using Plots # graph
using BenchmarkTools
# using BenchmarkProfiles
using MiniLogging

# basic configuration. The root logger level is then INFO
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")
loggerCR = get_logger("CR")
loggerCG = get_logger("CG")
# loggerminres = get_logger("minres")
loggerTR = get_logger("TrustRegion")
loggerTRLS = get_logger("TrustRegionLS")
loggerCR.level = MiniLogging.ERROR
loggerTR.level = MiniLogging.INFO
loggerCG.level = MiniLogging.ERROR
# loggerminres.level = MiniLogging.ERROR
loggerTRLS.level = MiniLogging.INFO

include("CG.jl")
include("CR.jl")
# include("minres.jl")
include("TrustRegion.jl")
include("TrustRegionLS.jl")
include("TRCR.jl")
include("TRCG.jl")
# include("TRminres.jl")
