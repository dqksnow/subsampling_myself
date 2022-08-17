function gendat(n, case, beta_t)
    d = length(beta_t)
    ds = d - 1
    corr  = 0.5
    sigmax = corr.^abs.(1(1:ds) .- (1:ds)')
    if case == 1 # Normal
        X = rand(MvNormal(ones(ds), sigmax), n)
    elseif case == 2 # lognormal
        X = exp.(rand(MvNormal(ones(ds), sigmax), n))
    elseif case == 3 # T3
        df = 5
        X = rand(MvNormal(ones(ds), sigmax), n)
        X = X ./ sqrt.(rand(Chisq(df), 1, n)./df)
    elseif case == 4 || case == 5 # Mix with order
        n5 = n ÷ 5
        X1 = rand(MvNormal(zeros(ds), sigmax), n5)
        X2 = rand(MvNormal(zeros(ds), sigmax), n5)
        X2 = X2 ./ sqrt.(rand(Chisq(2), 1, n5)./2)
        X3 = rand(MvNormal(zeros(ds), sigmax), n5)
        X3 = X3 ./ sqrt.(rand(Chisq(3), 1, n5)./3)
        X4 = 2rand(ds, n5)
        X5 = rand(MvNormal(zeros(ds), sigmax), n5)
        X5 = exp.(X5)
        X  = [X1 X2 X3 X4 X5]
        if case == 5 # Mix without order (random order)
            X = X[:, shuffle(1:end)]
        end
    elseif case == 6 # uniform
        X = randexp(ds, n) .- randexp(ds, n)

    end
    Y = beta_t[1] .+ X' * beta_t[2:end] + 3randn(n)
    Z = [ones(1, n); X]
    beta_f = (Z * Z') \ (Z * Y) 
    X̄ = mean(X, dims=2)
    Ȳ = mean(Y)

    F = svd(Z)
    # beta_f = F.V * Diagonal(1 ./ F.S) * F.U' * Y
    lv = sum(F.Vt.^2, dims=1)
    lv = lv ./ sum(lv)
    return X, Z, Y, X̄, Ȳ, beta_f, lv
end
