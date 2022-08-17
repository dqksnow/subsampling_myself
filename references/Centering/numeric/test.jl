initime = time()
using LinearAlgebra, Random, DelimitedFiles, Distributions

include("getidx.jl")
include("gendat.jl")

Random.seed!(1)
n = 10^5
d = 20
beta_t = -ones(d)
case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 2

# function simu(n, r, rpt, case, beta_t, nmd)
d = length(beta_t)
ds = d - 1
corr  = 0.5
sigmax = [corr+(1-corr)*(i==j) for i in 1:ds, j in 1:ds];
X, Z, Y, X̄, Ȳ, beta_f, lv = gendat(n, case, beta_t);

Xc = X .- X̄;
Yc = Y .- Ȳ;

# # opt
# e = Y .- Z'beta_f
# opt = sqrt.(lv) .* abs.(e')
# opt = opt ./ sum(opt)
# W_opt = 1 ./ opt; W_opt = W_opt / sum(W_opt)
# X̄_opt = sum(X ./ opt, dims=2) ./ sum(1 ./ opt)
# Ȳ_opt = sum(Y ./ opt') / sum(1 ./ opt)
# Xc_opt = X .- X̄_opt
# Yc_opt = Y .- Ȳ_opt

# A_opt = (Xc_opt * Xc_opt') \ ((Xc_opt .* W_opt) * Xc_opt') / (Xc_opt * Xc_opt')
# X̄_opt'A_opt*X̄_opt
# X̄'A_opt*X̄


# lev
W_lv = 1 ./ lv; W_lv = W_lv / sum(W_lv)
X̄_lv = sum(X ./ lv, dims=2) ./ sum(1 ./ lv)
Ȳ_lv = sum(Y ./ lv') / sum(1 ./ lv)
Xc_lv = X .- X̄_lv
Yc_lv = Y .- Ȳ_lv

A_lv = (Xc_lv * Xc_lv') \ ((Xc_lv .* W_lv) * Xc_lv') / (Xc_lv * Xc_lv')
X̄_lv'A_lv*X̄_lv
X̄'A_lv*X̄

B_lv = (Xc * Xc') \ ((Xc .* W_lv) * Xc') / (Xc * Xc')
(2/n) .* (X̄' / (Xc * Xc')) * (X̄ .- X̄_lv)

X̄_lv'B_lv*X̄_lv
X̄'B_lv*X̄

X̄_lv * X̄_lv'
X̄ * X̄'

X̄_lv'X̄_lv
X̄'X̄
