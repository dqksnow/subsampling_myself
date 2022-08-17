initime = time()
using DelimitedFiles, Random, Statistics, LinearAlgebra
using Distributions, Plots, LaTeXStrings
# using DataFrames 
# Threads.nthreads()
# using LinearAlgebra
BLAS.set_num_threads(1)
# ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
if "output" ∉ readdir() mkdir("output") end
if "csv" ∉ readdir("./output/") mkdir("./output/csv") end
include("gendat.jl")
include("estimators.jl")

Random.seed!(1)
N = 10^5
d = 10
beta0 = ones(d)
nmd = 5
case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 6
rpt = 1000
n = 1000
n0 = 1000
nss = [2000, 5000, 10000, 20000, 50000]
lns = length(nss)
Betas = fill(NaN, d, lns, rpt, nmd)

# (Z, Y) = gendat(N, case, beta0=beta0)
# X = [ones(N) Z]
# @time Thin(X, Y, nss, n0=1000);
# @time iboss(Z, Y, nss);
# @time Uni(X, Y, nss.+n0, poi=false);
# @time Poi(X, Y, nss, n0=1000, estH=true);
# @time Rep(X, Y, nss, n0=1000);

@time for rr in 1:rpt
    (Z, Y) = gendat(N, case, beta0=beta0)
    X = [ones(N) Z]
    m = 0
    Betas[:,:,rr,m+=1] = Thin(X, Y, nss)
    Betas[:,:,rr,m+=1] = Uni(X, Y, nss)
    Betas[:,:,rr,m+=1] = iboss(Z, Y, nss)
    Betas[:,:,rr,m+=1] = Poi(X, Y, nss, n0=n0, estH=true)
    Betas[:,:,rr,m+=1] = Rep(X, Y, nss, n0=n0)
end

print("case:$(case)\n N:$(N)\n rpt:$(rpt) \n nss:$(nss)\n")
rec = fill(NaN, nmd, length(nss))
# @time cal = simu(X, Y, nss, rpt, nmd)
for m in 1:nmd, (idn, n) in enumerate(nss)
    rec[m,idn] = sum(mean((Betas[:,idn,:,m] .- beta0).^2, dims=2))
    fname = "output/csv/case$(case)N$(N)method$(m)n$(n).csv"
    writedlm(fname, Betas[:,idn,:,m])
end

show(stdout, "text/plain", rec)
print("\n", sum(rec))
# label = ["dc: M=" .* map(string, nss); "uniform"]
label = ["Thin" "R uni" "iboss" "P opt h" "P uni"]
pl = plot(nss./N, log.(rec'), label=label, lw=2, m=(7,:auto))
# pl = plot(nss,   log10.(rec'), label=label, lw=2, m=(7,:auto))
# pl = plot(log10.(nss), log10.(rec'), label=label, lw=2, m=(7,:auto))
# plot(log10.(nss), log10.(rec'), label=label, lw=2, m=(7,:auto))
savefig(pl, "output/0case$(case).pdf")
writedlm("output/case$(case).csv", [nss rec'])

print("\n Total time used: $(round(time() - initime)) seconds")

# scatter!(nsss, rec', legend=false)
# scatter!(nsss, rec', legend=false)

# rec = Array{Float64}(undef, 3, 2)
# # mapslices
# using Statistics
# nanmean(x) = mean(filter(!isnan,x))
# nanmean(x,y) = mapslices(nanmean,x,dims=y)
# y = [NaN 2 3 4;5 6 NaN 8;9 10 11 12]
# nanmean(y)
# nanmean(y,1)
