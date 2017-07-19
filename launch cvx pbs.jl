Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

for problem in Problems
  model = MathProgNLPModel(problem(), name=string(problem))
  TRCG(model)
  reset!(model)
  TRCR(model)
  reset!(model)
  TRminres(model)
end
