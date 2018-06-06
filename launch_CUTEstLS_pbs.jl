using NLPModels
using CUTEst
using Krylov
include("TRLS.jl")

num = parse(Int64,ARGS[1])
Problems = open(readlines, "ProblemsLS.txt")

Algos = [lsqr, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #r  #Av  #A'v  #it #s_it #vs_it #rej_it optimal"]
# @info(loggerCUTEstLS, D[1])

problem = Problems[num]
model = CUTEstModel(problem)
lsmodel = FeasibilityResidual(model)
for algo in Algos
    T = TRLS(lsmodel, algo)
    S = @sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d %5s", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14], T[15])
    # @info(loggerCUTEstLS, S)
    D = vcat(D, S)
    reset!(lsmodel)
end
finalize(model)
writedlm(string(problem, "_CUTEstLS.txt"), D)
