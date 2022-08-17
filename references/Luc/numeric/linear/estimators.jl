# initime = time()
# using Distributions, Random, Statistics, Plots, LaTeXStrings, LinearAlgebra
# using DelimitedFiles
# # using DataFrames 
# # Threads.nthreads()
# # using LinearAlgebra
# # BLAS.set_num_threads(1)
# # ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
# include("gendat.jl")
# include("seqD.jl")
# include("getidx.jl")

# Random.seed!(0)
# case = 1
# N = 10^5
# (Z, Y) = gendat(N, case, beta0=ones(10))
# X = [ones(N) Z]

function seqD(X::Matrix{Float64}, n::Int=10^4;
                k0::Int=5size(X)[2], q::Float64=0.625, gamma::Float64=0.1)
    N = size(X)[1]
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
    return idx, nk
end
# idx = seqD(X, n)


function Thin(X::Matrix{Float64}, Y::Vector{Float64}, nss::Vector{Int};
                n0::Int=1000, q::Float64=0.625, gamma::Float64=0.1)
    (N,d) = size(X)
    lns = length(nss)
    Betas = fill(NaN, d, lns)    
    for (idn, n) in enumerate(nss)
        # print("n=$(n)\n")
        idx_thin = seqD(X, n, k0=n0, q=0.625, gamma=0.1)[1]
        x_thin = X[idx_thin, :]
        y_thin = Y[idx_thin]
        beta_thin = (x_thin'x_thin) \ x_thin'*y_thin
        Betas[:,idn] = beta_thin
    end
    return Betas
end
# @time sum(Thin(X, Y, nss, n0=1000))

function getidx(x::Matrix, k::Int64)
    (n,p) = size(x)
    idx = Array{Int64}(undef, k)
    counter = 1
    r = cld(k, 2p)
    tmp2 = x[:,1]
    l = sort!(tmp2; alg=PartialQuickSort(r))[r]
    u = sort!(tmp2; alg=PartialQuickSort(r), rev=true)[r]
    bl = 1; bu = 1;
    for s in 1:n
        # global counter
        if bl > r && bu > r
            break
        elseif bl <= r && x[s,1] <= l
            idx[counter] = s
            counter += 1
            bl += 1
        elseif bu <= r && x[s,1] >= u
            idx[counter] = s
            counter += 1
            bu += 1
        end
    end
    for j in 2:p
        tl = n - counter + 1
        tmp = Vector{Float64}(undef, tl)
        t = 1; v = 1;
        for s in 1:n
            # global v, t
            if s != idx[v]
                tmp[t] = x[s,j]
                t += 1
            elseif v < counter - 1
                v += 1
            end
        end
        l = sort!(tmp; alg=PartialQuickSort(r))[r]
        u = sort!(tmp; alg=PartialQuickSort(r), rev=true)[r]
        v = 1; bl = 1; bu = 1;
        for s in 1:n
            # global v, k, counter
            if counter > k || (bl > r && bu > r)
                break
            elseif v < counter && s == idx[v]
                v += 1
            elseif bl <= r && x[s,j] <= l
                idx[v+1:counter] = idx[v:counter-1]
                idx[v] = s
                v += 1
                counter += 1
                bl += 1
            elseif bu <= r && x[s,j] >= u
                idx[v+1:counter] = idx[v:counter-1]
                idx[v] = s
                v += 1
                counter += 1
                bu += 1
            end
        end
    end
    return idx
end
# @time new = sum(getidx(Z, 10000));

# IBOSS
function iboss(Z::Matrix{Float64}, Y::Vector{Float64}, nss::Vector{Int})
    (N,d) = size(Z) .+ (0,1)
    lns = length(nss)
    Betas = fill(NaN, d, lns)    
    for (idn, n) in enumerate(nss)
        # print("n=$(n)\n")
        idx_iboss = getidx(Z, n)
        x_iboss = [ones(n) Z[idx_iboss, :]]
        y_iboss = Y[idx_iboss]
        beta_iboss = (x_iboss'x_iboss) \ x_iboss'*y_iboss
        Betas[:,idn] = beta_iboss
    end
    return Betas
end
# iboss(Z, Y, nss)


# Poisson opt
function Poi(X::Matrix{Float64}, Y::Vector{Float64}, nss::Vector{Int};
                n0::Int=1000, a::Float64=0.05, estH=false)
    (N,d) = size(X)
    lns = length(nss)
    Betas = fill(NaN, d, lns)    
    # idx_plt = sample(1:N, n0, replace=false)
    u_plt = rand(N)
    idx_plt = u_plt .<= n0/N 
    x_plt = X[idx_plt, :]
    y_plt = Y[idx_plt]
    n0s = length(y_plt)
    ddm_plt = (x_plt'x_plt) 
    beta_plt = ddm_plt \ x_plt'*y_plt
    e_plt = Y .- X * beta_plt
    dm = abs.(e_plt) .* sqrt.(vec(sum(X.^2, dims=2)))
    pi_P = (1-a) .* dm ./ (N*mean(dm[idx_plt])) .+ a/N
    if estH
        dmh = min.(dm, quantile(dm, 1-n/3N))
        pi_P = (1-a) .* dmh ./ (N*mean(dmh[idx_plt])*(n0s/(n0s-d))) .+ a/N
    end
    for (idn, n) in enumerate(nss)
        # print("n=$(n)\n")
        u = rand(N)
        idx_popt = u .<= (n-n0).*pi_P
        x_popt = X[idx_popt, :]
        y_popt = Y[idx_popt]
        pi_popt = min.((n-n0) .* pi_P[idx_popt], 1)
        ddm_popt = (x_popt' * (x_popt ./ pi_popt))
        beta_popt = ddm_popt \ (x_popt ./ pi_popt)'*y_popt
        beta_popt = (ddm_popt ./(N/(n-n0)) .+ ddm_plt) \ (ddm_popt*beta_popt ./(N/(n-n0)) .+ ddm_plt*beta_plt)
        Betas[:,idn] = beta_popt
    end
    return Betas
end
# Poi(X, Y, nss, n0=1000, estH=true)

# Sampling with replacement 
function Rep(X::Matrix{Float64}, Y::Vector{Float64}, nss::Vector{Int};
                n0::Int=1000, a::Float64=0.05)
    (N,d) = size(X)
    lns = length(nss)
    Betas = fill(NaN, d, lns)
    idx_plt = sample(1:N, n0, replace=false)
    x_plt = X[idx_plt, :]
    y_plt = Y[idx_plt]
    ddm_plt = (x_plt'x_plt) 
    beta_plt = ddm_plt \ x_plt'*y_plt
    e_plt = Y .- X * beta_plt
    dm = abs.(e_plt) .* sqrt.(vec(sum(X.^2, dims=2)))
    pi_R = (1-a) .* dm ./ sum(dm) .+ a/N
    for (idn, n) in enumerate(nss)
        # print("n=$(n)\n")
        idx_ropt = wsample(1:N, pi_R, n-n0, replace=true)
        x_ropt = X[idx_ropt, :]
        y_ropt = Y[idx_ropt]
        pi_ropt = pi_R[idx_ropt]
        ddm_ropt = (x_ropt' * (x_ropt ./ pi_ropt))
        beta_ropt = ddm_ropt \ (x_ropt ./ pi_ropt)'*y_ropt
        beta_ropt = (ddm_ropt ./N .+ ddm_plt) \ (ddm_ropt*beta_ropt ./N .+ ddm_plt*beta_plt)
        Betas[:,idn] = beta_ropt
    end
    return Betas
end
# Rep(X, Y, nss, n0=1000)

# Uniform sampling
function Uni(X::Matrix{Float64}, Y::Vector{Float64}, nss::Vector{Int};
                poi=true)
    (N,d) = size(X)
    lns = length(nss)
    Betas = fill(NaN, d, lns)
    for (idn, n) in enumerate(nss)
        # print("n=$(n)\n")
        if poi
            u = rand(N)
            idx_uni = u .<= n/N
        else
            idx_uni = sample(1:N, n, replace=true)
        end
        x_uni = X[idx_uni, :]
        y_uni = Y[idx_uni]
        beta_uni = (x_uni'x_uni) \ x_uni'*y_uni
        Betas[:,idn] = beta_uni
    end
    return Betas
end
# Uni(X, Y, nss, poi=false)

# Random.seed!(1)
# N = 10^5
# d = 51
# beta0 = ones(d)
# nmd = 5
# case = haskey(ENV, "case") ? parse(Int, ENV["case"]) : 1
# rpt = 10
# n = 1000
# nss = [2000, 5000, 10000, 20000, 50000]

# (X, Y, beta_f, lv) = gendat(N, case, beta0)

# print("case:$(case)\n N:$(N)\n rpt:$(rpt) \n nss:$(nss)\n")
# rec = fill(NaN, nmd, length(nss))
# @time cal = simu(X, Y, nss, rpt, nmd)
# for m in 1:nmd, idn in 1:length(nss)
#     rec[m,idn] = sum(mean((cal[:,:,m, idn] .- beta_f).^2, dims=2))
#     fname = "output/case$(case)N$(N)method$(m)n$(n).csv"
#     writedlm(fname, cal[:,:,m,idn])
# end

# show(stdout, "text/plain", rec)
# print("\n", sum(rec))
# # label = ["dc: M=" .* map(string, nss); "uniform"]
# label = ["R opt", "R uni", "P opt", "P opt h", "P uni"]
# pl = plot(nss/N, log.(rec'), label=label, lw=2, m=(7,:auto))
# plot(nss, log.(rec'), label=label, lw=2, m=(7,:auto))
# # pl = plot(log10.(nss), log10.(rec'), label=label, lw=2, m=(7,:auto))
# # plot(log10.(nss), log10.(rec'), label=label, lw=2, m=(7,:auto))
# savefig(pl, "output/case$(case).pdf")
# writedlm("output/case$(case).csv", [nss rec'])

# print("\n Total time used: $(round(time() - initime)) seconds")

# # scatter!(nss, rec', legend=false)
# # scatter!(nss, rec', legend=false)

# # rec = Array{Float64}(undef, 3, 2)
# # # mapslices
# # using Statistics
# # nanmean(x) = mean(filter(!isnan,x))
# # nanmean(x,y) = mapslices(nanmean,x,dims=y)
# # y = [NaN 2 3 4;5 6 NaN 8;9 10 11 12]
# # nanmean(y)
# # nanmean(y,1)

function seqPhi(X::Matrix{Float64}, n::Int=10^4;
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
    Phi = Vector{Float64}(undef, N)
    Phi[1:k0] .= log(det(M0))
    nk = k0; Mk = M0; Ck = C0; fk = f0;
    for k in (k0+1):N
        # global nk, Mk, Ck, fk
        # push!(Phi, log(det(Mk)))
        Phi[k] = log(det(Mk))
        xkn = X[k,:]
        if nk >= n
            Phi[idx[nk]+1:N] .= log(det(Mk))
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
    return idx, nk, k0, Phi
end
