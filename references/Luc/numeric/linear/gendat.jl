# using Feather

function gendat(n::Int, case::Int=1; beta0::Vector=ones(10))
    d = length(beta0)
    ds = d - 1
    corr  = 0.5
    sigmax = [corr+(1-corr)*(i==j) for i in 1:ds, j in 1:ds]
    if case == 1 # Normal
        Z = collect(rand(MvNormal(zeros(ds), sigmax), n)')
    elseif case == 2 # lognormal
        Z = exp.(rand(MvNormal(zeros(ds), sigmax), n)')
    elseif case == 3 # T2
        df = 2
        Z = rand(MvNormal(zeros(ds), sigmax), n);
        Z = collect(Z') ./ sqrt.(rand(Chisq(df), n)./df)
    elseif case == 4 # T3
        df = 3
        Z = rand(MvNormal(zeros(ds), sigmax), n);
        Z = collect(Z') ./ sqrt.(rand(Chisq(df), n)./df)
    # elseif case == 4 # Mix with order
    #     n5 = cld(n, 5)
    #     Z1 = rand(MvNormal(zeros(ds), sigmax), n5)'
    #     Z2 = rand(MvNormal(zeros(ds), sigmax), n5)
    #     Z2 = collect(Z2') ./ sqrt.(rand(Chisq(2), n5)./2)
    #     Z3 = rand(MvNormal(zeros(ds), sigmax), n5)
    #     Z3 = collect(Z3') ./ sqrt.(rand(Chisq(3), n5)./3)
    #     Z4 = 2rand(n5, ds)
    #     Z5 = rand(MvNormal(zeros(ds), sigmax), n5)'
    #     Z5 = exp.(Z5)
    #     Z  = [Z1; Z2; Z3; Z4; Z5]
    #     # mix = MixtureModel(MvNormal,
    #     #     [(ones(ds), sigmax), (-ones(ds), sigmax)],
    #     #     [0.5, 0.5])
    #     # Z = rand(mix, n);
    #     # Z = convert(Array{Float64,2}, Z')
    elseif case == 5 # Mix without order (random order)
        n5 = cld(n, 5)
        Z1 = rand(MvNormal(zeros(ds), sigmax), n5)'
        Z2 = rand(MvNormal(zeros(ds), sigmax), n5)
        Z2 = collect(Z2') ./ sqrt.(rand(Chisq(2), n5)./2)
        Z3 = rand(MvNormal(zeros(ds), sigmax), n5)
        Z3 = collect(Z3') ./ sqrt.(rand(Chisq(3), n5)./3)
        Z4 = 2rand(n5, ds)
        Z5 = rand(MvNormal(zeros(ds), sigmax), n5)'
        Z5 = exp.(Z5)
        Z  = [Z1; Z2; Z3; Z4; Z5]
        Z = Z[shuffle(1:end), :]
    elseif case == 6 # exponential
        Z = randexp(n, ds)
        # elseif case == 6 # uniform
        #     Z = rand(n, ds)
    # elseif case == 10 # sensor data
    #     # dat = readdlm("/home/ossifragus/ondisk/ethylene_CO.txt")
    #     datf = Feather.read("/home/ossifragus/ondisk/ethylene_CO.feather")
    #     dat = convert(Array{Float64,2}, datf)
    #     Y = log.(dat[20001:end, 19])
    #     Z = log.(dat[20001:end, [4; 6:18]])
    #     n = length(Y)
    end
    if case != 10
        Y = beta0[1] .+ Z * beta0[2:end] + randn(n)
    end
    # X = [ones(n) Z]

    # F = svd(X)
    # beta_f = F.V * Diagonal(1 ./ F.S) * F.U' * Y
    # lv = sqrt.(vec(sum(X.^2, dims=2)))
    # # e = Y - X * beta_f
    # # lv = abs.(e) .* sqrt.(vec(sum(X.^2, dims=2)))
    # # lv = vec(sum(F.U.^2, dims=2))
    # # lv = lv ./ sum(lv)

    return Z, Y# , beta_f, lv
end
