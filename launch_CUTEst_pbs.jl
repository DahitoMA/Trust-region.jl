using CUTEst
using Optimize

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
Problems = deleteat!(Problems,56) # pb "INDEF" not treated
Algos = [CG, CR]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it"]

for problem in Problems
    model = CUTEstModel(problem)
    for algo in Algos
        # (fx, fx0, gNorm, g0Norm, iter, status, nf, ng, nhv, s_iter, vs_iter, r_iter) = TrustRegion(model, algo)
        # @printf("%-15s  %5s %8d  %9.2e  %9.2e  %7.1e  %7.1e  %5d  %5d  %5d  %6d  %6d  %6d  %6d  %s\n",
        #     model.meta.name, string(algo), model.meta.nvar, fx, fx0, gNorm, g0Norm,
        #     nf, ng, nhv, iter, s_iter, vs_iter, r_iter, status)

        # T = TrustRegion(model, algo)
        D = vcat(D, TrustRegion(model, algo))
        reset!(model)
    end
    finalize(model)
end

writedlm("CUTEstpbs.txt", D)
