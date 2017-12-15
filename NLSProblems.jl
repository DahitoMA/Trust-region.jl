using NLSProblems
using Krylov

Problems = [mgh01, mgh02, mgh03,mgh04, mgh05, mgh06, mgh07, mgh08, mgh09,
mgh10, mgh11, mgh12, mgh13, mgh14, mgh15, mgh16, mgh17, mgh18, mgh19,
mgh20, mgh21, mgh22, mgh23, mgh24, mgh25, mgh26, mgh27, mgh28, mgh29,
mgh30, mgh31, mgh32, mgh33, mgh34, mgh35]

Algos = [cgls, lsqr, crls, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #r  #Av  #A'v  #it #s_it #vs_it #rej_it"]

for problem in Problems
    model = ADNLSModel(problem())
    for algo in Algos
        T = TrustRegionLS(model, algo)
        D = vcat(D, T)
        println(T)
        reset!(model)
    end
end

writedlm("NLSProblems.txt", D)

lines = open(readlines, "NLSProblems.txt")
V = []
i = 0
for l in lines
    i += 1
    l = replace(l, "Generic", string("mgh", Int(ceil((i-1)/4))))
    V = vcat(V, l)
end
writedlm("NLSProblems.txt", V)
