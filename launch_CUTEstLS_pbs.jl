using CUTEst
using Krylov

Problems = open(readlines, "ProblemsLS.txt")
Problems = Problems[1:22]
Algos = [cgls]#, lsqr, crls, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it"]

for problem in Problems
    model = CUTEstModel(problem)
    for algo in Algos
        D = vcat(D, TrustRegionLS(model, algo))
        println(D)
        reset!(model)
    end
    finalize(model)
end

writedlm("CUTEstLSpbs.txt", D)
