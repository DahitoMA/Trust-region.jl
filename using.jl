using Lint # check code
using JLD # save variables
using OptimizationProblems
using NLPModels
using LinearOperators
using PyPlot # graph
using BenchmarkTools
using BenchmarkProfiles

include("CG.jl")
include("CR.jl")
include("minres.jl")
include("TRCG.jl")
include("TRCR.jl")
include("TRminres.jl")
include("TrustRegion.jl")
include("CRtrunc.jl")

using Krylov
