initime = time()
using Distributions, Random, Statistics, Plots, LaTeXStrings, LinearAlgebra
using DelimitedFiles
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
case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 1
rpt = 1000
n = 10000
n0 = 1000
nss = [2000, 5000, 10000, 20000, 50000]
lns = length(nss)
Betas = fill(NaN, d, lns, rpt, nmd)

(Z, Y) = gendat(N, case, beta0=beta0)
# Z = randn(N); Z = [Z Z.^2];
X = [ones(N) Z]

function seqD(X::Matrix{Float64}, n::Int=10^4;
                k0::Int=5size(X)[2], q::Float64=0.625, gamma::Float64=0.1)
    (N, d) = size(X)
    a = n/N
    x0 = X[1:k0,:]
    M0 = x0'x0 ./ k0
    M0iv = inv(M0)
    zt0 = Vector{Float64}(undef, k0)
    for i in 1:k0
        zt0[i] = x0[i,:]' * M0iv * x0[i,:] # tr(M0iv * (x0[i,:]x0[i,:]'))
    end
    sort!(zt0)
    k0p = ceil(Int, (1-a/2)*k0)
    k0m = floor(Int, (1-3a/2)*k0)
    C0 = zt0[ceil(Int, (1-a)*k0)]
    b0 = k0 / (k0p-k0m)
    h = zt0[k0p] - zt0[k0m]
    h0 = h / k0^gamma
    f0 = sum(abs.(zt0 .- C0) .<= h0) / 2k0*h0
    idx = Vector{Int}(undef, n)
    idx[1:k0] = 1:k0
    # s = Array{Float64}(undef, n, d)
    # s[1:k0,:] = x0
    # Phi = Float64[]
    nk = k0; Mk = M0; Ck = C0; fk = f0;
    for k in (k0+1):N
        # global nk, Mk, Ck, fk
        # push!(Phi, log(det(Mk)))
        # xkn = [1; X[k,:]]
        xkn = X[k,:]
        if nk >= n
            break
        elseif n-nk >= N-k
            zkp = Inf
        else
            zkp = xkn' / Mk * xkn # tr(inv(Mk) * (xkn*xkn'))
        end
        if zkp > Ck
            nk += 1
            Mk = (nk-1)/nk .* Mk + xkn*xkn' ./nk
            # s[nk,:] = [1; X[k,:]]
            idx[nk] = k
        end
        bk = min(1/fk, b0*k^gamma)
        Ck = Ck + bk/(k+1)^q * ((zkp > Ck)-a)
        hkp = h / k^gamma
        fk = fk + 1/k^q * ((abs(zkp-Ck)<=hkp)/(2hkp) - fk)
    end
    return idx, nk, k0
end
(idx, nk, k0) = seqD(X, n)# , k0=n0)

# s = [ones(n) X[idx[1],:]]

for (icase, case) in enumerate(1:5)
(Z, Y) = gendat(N, case, beta0=beta0)
# Z = randn(N); Z = [Z Z.^2];
X = [ones(N) Z]

# Calculate Phi
x0 = X[1:k0,:]
Mk = x0'x0 ./ k0
Phi = Vector{Float64}(undef, N)
Phi[1:k0] .= log(det(Mk))
nk = k0+1
for k in (k0+1):N
    global Mk, nk
    Phi[k] = log(det(Mk))
    # push!(Phi, log(det(Mk)))
    if nk >n
        Phi[idx[nk-1]+1:N] .= log(det(Mk))
        break
    elseif k == idx[nk]
        xkn = X[k,:]
        Mk = (nk-1)/nk .* Mk + xkn*xkn' ./nk
        nk += 1
    end
end

plot(log10.(k0:length(Phi)), Phi[k0:N], legend=false, xlim=[1, log10.(N)],
     xlabel=L"\mathbf{\log_{10}(k)}", xguidefontsize=15,
     ylabel=L"\mathbf{\Phi(M_{n_k})}", yguidefontsize=15)
# hline!([1.6354, 3.2963])
