initime = time()
using DelimitedFiles, Random, Statistics, LinearAlgebra
using Distributions, Plots, LaTeXStrings
# using DataFrames 
# Threads.nthreads()
# using LinearAlgebra
BLAS.set_num_threads(1)
# ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
if "output" ∉ readdir() mkdir("output") end
if "phi" ∉ readdir("./output/") mkdir("./output/phi") end
include("gendat.jl")
include("estimators.jl")

Random.seed!(1)
N = 10^5
d = 10
beta0 = ones(d)
nmd = 5
case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 6
rpt = 10
n = 2000
n0 = 1000
k0 = 5d
nss = [2000, 5000, 10000, 20000, 50000]
lns = length(nss)
Betas = fill(NaN, d, lns, rpt, nmd)

label = reshape("n=" .* map(string, nss), 1, nmd)
case_ss = collect(1:6)
stp = ["pdf", "png"][2]
for i in 1:length(case_ss)
    case = case_ss[i]
    PhiT = Array{Float64}(undef, N, lns, rpt)
    for j in 1:rpt
        (Z, Y) = gendat(N, case, beta0=beta0)
        X = [ones(N) Z]
        PhiM = Array{Float64}(undef, N, lns)
        for (idn, n) in enumerate(nss)
            PhiT[:,idn,j] = seqPhi(X, n)[4]
            # (idx, nk, k0, Phi) = seqPhi(X, n)
            # PhiM[:,idn] = Phi
        end
    end
    PhiM = mean(PhiT, dims=3)
    pl = plot(log10.(k0:N), PhiM[k0:N,:], legend=:topleft,
              xlim=[log10.(k0), log10.(N)], label=label,
              xlabel=L"\mathbf{\log_{10}(k)}", xguidefontsize=15,
              ylabel=L"\mathbf{\Phi(M_{n_k})}", yguidefontsize=15)
    # hline!([1.6354, 3.2963])
    # savefig(pl, "output/phi/PhiCase$(case).pdf")
    savefig(pl, "output/phi/PhiCase$(case).$(stp)")
end
print("\n Total time used: $(round(time() - initime)) seconds")

cbfile = "./output/0Phi_linear.pdf"
fns = "output/phi/PhiCase" .* map(string, case_ss) .* ".$(stp)"
stp == "pdf" ? run(`pdftk $(fns) output $(cbfile)`) :
    run(`convert $(fns)  $(cbfile)`)
run(`cp $(cbfile) ../draft/figures/`)
