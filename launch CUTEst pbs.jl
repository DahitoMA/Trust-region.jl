using CUTEst

Problems = open(readlines, "cutest_unc.lst")

# pb_to_skip = []
# for pb in Problems
#     nlp = CUTEstModel(pb)
#     if nlp.meta.nvar < 10
#         push!(pb_to_skip, pb)
#     end
#     CUTEst.cutest_finalize(nlp)
# end

pb_to_skip = open(readlines, "pb_to_skip.txt")
Problems = filter(x -> !in(x, pb_to_skip), Problems)

Algos = [CG, CR]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it"]

for problem in Problems
    for algo in Algos
        model = CUTEstModel(problem)
        D = vcat(D, TrustRegion(model, algo))
        CUTEst.cutest_finalize(model)
    end
end

writedlm("CUTEstpbs.txt", D)
