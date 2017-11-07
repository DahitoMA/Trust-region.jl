Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

Algos = [CG, CR, minres]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it"]

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    for algo in Algos
        D = vcat(D, TrustRegion(model, algo))
        reset!(model)
    end
end

# writedlm("cvxpbs2.txt", D)
writedlm("cvxpbsmincg.txt", D)
