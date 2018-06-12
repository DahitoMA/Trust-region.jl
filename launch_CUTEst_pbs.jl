using NLPModels
using CUTEst
include("TrustRegionCUTEst.jl")
include("CG.jl")
include("CR.jl")

# using Optimize

# Problems = open(readlines, "cutest_unc.lst")

# pb_to_skip = []
# for pb in Problems
#     nlp = CUTEstModel(pb)
#     if nlp.meta.nvar < 10
#         push!(pb_to_skip, pb)
#     end
#     CUTEst.cutest_finalize(nlp)
# end

# pb_to_skip = open(readlines, "pb_to_skip.txt")
# Problems = filter(x -> !in(x, pb_to_skip), Problems)
# writedlm("Problems.txt", Problems)

num = parse(Int64,ARGS[1])
Problems = open(readlines, "Problems.txt")
Algos = [CG, CR]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it optimal"]
# @info(loggerCUTEst, D[1])

problem = Problems[num]
model = CUTEstModel(problem)
for algo in Algos
    T = TrustRegionCUTEst(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %4d %4d %4d %5s", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14], T[15])
    # @info(loggerCUTEst, S)
    D = vcat(D, S)
    reset!(model)
end
finalize(model)
writedlm(string(problem, "_CUTEst.txt"), D)
