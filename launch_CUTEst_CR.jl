using CUTEst
import CSV
import DataFrames
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

Problems = open(readlines, "Problems.txt")
# Problems = deleteat!(Problems,56) # pb "INDEF" not treated
algo = CR
csvname = string("CUTEstTR",algo,".csv")

D = ["model" "nvar"   "f(x)"    "f(x0)"   "‖g(x)‖"  "‖g(x0)‖"   "#f"  "#g"  "#Hv"  "#it" "#s_it" "#vs_it" "#rej_it"]
@info(loggerCUTEstCR, D[1])

for problem in Problems
    model = CUTEstModel(problem)
    T = TrustRegionCUTEst(model, algo)
    S = @sprintf("%5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %4d %4d %4d", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13])
    @info(loggerCUTEstCR, S)
    D = vcat(D, T)
    reset!(model)
    finalize(model)
end

DF = DataFrames.DataFrame(D)
CSV.write(csvname, DF)
