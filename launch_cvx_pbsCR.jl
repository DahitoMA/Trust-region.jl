using OptimizationProblems
import CSV
import DataFrames

Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

algo = CR
csvname = string("cvxTR",algo,".csv")

D = ["model" "nvar"   "f(x)"    "f(x0)"   "‖g(x)‖"  "‖g(x0)‖"   "#f"  "#g"  "#Hv"  "#it" "#s_it" "#vs_it" "#rej_it"]
@info(loggercvx, D[1])

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    T = TrustRegion(model, algo)
    S = @sprintf("%5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %4d %4d %4d", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13])
    @info(loggercvx, S)
    D = vcat(D, T)
    reset!(model)
end
DF = DataFrames.DataFrame(D)
CSV.write(csvname, DF)
