using OptimizationProblems

Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

# Algos = [CG, CR, minres]
Algos = [CG, CR]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it"]
@info(loggercvx, D[1])

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    for algo in Algos
        T = TrustRegion(model, algo)
        S = @sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14])
        @info(loggercvx, S)
        D = vcat(D, S)
        reset!(model)
    end
end

writedlm("cvxpbsCGCR.txt", D)
