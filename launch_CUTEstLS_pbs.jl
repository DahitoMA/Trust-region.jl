using CUTEst
using Krylov

Problems = open(readlines, "ProblemsLS.txt")

Algos = [cgls, lsqr, crls, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #r  #Av  #A'v  #it #s_it #vs_it #rej_it"]
@info(loggerCUTEstLS, D[1])

for problem in Problems
    model = CUTEstModel(problem)
    lsmodel = FeasibilityResidual(model)
    for algo in Algos
        T = TrustRegionLS(lsmodel, algo)
        S = @sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14])
        @info(loggerCUTEstLS, S)
        D = vcat(D, S)
        reset!(lsmodel)
    end
    finalize(model)
end

writedlm("CUTEstLSpbs.txt", D)
