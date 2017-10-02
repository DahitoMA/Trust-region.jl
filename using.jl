# using Lint # check code
# using JLD # save variables
using OptimizationProblems
using NLPModels
# using LinearOperators
# using PyPlot # graph
# using BenchmarkTools
# using BenchmarkProfiles
using MiniLogging
# using Krylov

# basic configuration. The root logger level is then INFO
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")
logger = get_logger("CR") # get root logger


# include("CG.jl")
include("CR.jl")
# include("minres.jl")
include("TrustRegion.jl")
