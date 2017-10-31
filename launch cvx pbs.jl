Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

Algos = [CG, CR]#, minres]

loggerCvx = MiniLogging.Logger("loggerCvx", MiniLogging.INFO)
push!(loggerCvx.handlers, MiniLogging.Handler(open("cvxpbs.txt", "w"), "%Y-%m-%d %H:%M:%S"))

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    for algo in Algos
        @info(loggerCvx, TrustRegion(model, algo))
        reset!(model)
    end
end
