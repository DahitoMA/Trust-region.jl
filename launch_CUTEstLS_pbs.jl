using CUTEst
using Krylov

Problems = open(readlines, "ProblemsLS.txt")
Problems = Problems[8:9]
Algos = [cgls]#[cgls, lsqr, crls, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #r  #Av  #A'v  #it #s_it #vs_it #rej_it"]

for problem in Problems
    model = CUTEstModel(problem)
    lsmodel = FeasibilityResidual(model)
    for algo in Algos
        T = TrustRegionLS(lsmodel, algo)
        D = vcat(D, T)
        println(T)
        reset!(lsmodel)
    end
    finalize(model)
end

writedlm("CUTEstLSpbs.txt", D)
