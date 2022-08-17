initime = time()
using LinearAlgebra, Random, DelimitedFiles
using Distributions, ProgressMeter, FLoops
# using Plots, LaTeXStrings
BLAS.set_num_threads(1) # BLAS.get_num_threads()

if "output" ∉ readdir() mkdir("output") end
if "csv" ∉ readdir("./output/") mkdir("./output/csv") end

include("getidx.jl")
include("gendat.jl")
include("simu.jl")

Random.seed!(1)
n = 10^5
d = haskey(ENV, "dims") ? parse(Int, ENV["dims"]) : 20
beta_t = ones(d)
nmd = 9
case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 2
rpt = 1000
r = 100
print("d:$(d); case:$(case); rpt:$(rpt); r:$(r)\n")
nss = [5*10^3, 10^4, 10^5, 10^6][[3]]
rec = fill(NaN, nmd, length(nss))
rec0 = copy(rec); rec1 = copy(rec);

for (id, n) in enumerate(nss)
    print("n=$n: ")
    @time cal = simu(n, r, rpt, case, beta_t, nmd)
    for m in 1:nmd
        rec[m,id] = sum(mean((cal[:,:,m] .- beta_t).^2, dims=2)[1:end])
        rec0[m,id] = sum(mean((cal[:,:,m] .- beta_t).^2, dims=2)[1])
        rec1[m,id] = sum(mean((cal[:,:,m] .- beta_t).^2, dims=2)[2:end])
        fname = "output/csv/d$(d)case$(case)n$(n)method$(m).csv"
        writedlm(fname, cal[:,:,m])
    end
end

label = ["uni"  "uni-c" "iboss" "iboss-c" "lev" "lev-c" "opt" "opt-c" "full"]
# show(stdout, "text/plain", [vec(label) rec]); print("\n")
show(stdout, "text/plain", [vec(label) rec0]); print("\n")
show(stdout, "text/plain", [vec(label) rec1./(d-1)]); print("\n")
# show(stdout, "text/plain", [vec(label) rec][1:end,:])
print("\n", sum(rec))
# pl = plot(log10.(nss), log10.(rec'), label=label, lw=2, m=(7,:auto));
# savefig(pl, "output/case$(case).pdf")

h = ["Case$(case)" repeat(["Intercept" "Slope"], outer=(1, length(nss)))]
writedlm("output/d$(d)case$(case).csv", [h; vec(label) rec0 rec1./(d-1)])

print("\n Total time used: $(round(time() - initime)) seconds")

# scatter!(rss, rec', legend=false)
# scatter!(rss, rec', legend=false)

# rec = Array{Float64}(undef, 3, 2)
# # mapslices
# using Statistics
# nanmean(x) = mean(filter(!isnan,x))
# nanmean(x,y) = mapslices(nanmean,x,dims=y)
# y = [NaN 2 3 4;5 6 NaN 8;9 10 11 12]
# nanmean(y)
# nanmean(y,1)
