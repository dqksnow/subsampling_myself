function simu(n, r, rpt, case, beta_t, nmd)
    d = length(beta_t)
    ds = d - 1
    corr  = 0.5
    sigmax = [corr+(1-corr)*(i==j) for i in 1:ds, j in 1:ds]
    Betas = fill(NaN, d, rpt, nmd)
    # pgs = Progress(rpt)
    @floop ThreadedEx(basesize=rpt÷15) for rr in 1:rpt
    # Threads.@threads :static for rr in 1:rpt
    # for rr in 1:rpt
        X, Z, Y, X̄, Ȳ, beta_f, lv = gendat(n, case, beta_t)

        # uniform
        idx_uni = (1:n)[rand(n) .<= r/n]
        z_uni = Z[:, idx_uni]
        y_uni = Y[idx_uni]
        beta_uni = (z_uni * z_uni') \ (z_uni * y_uni)

        x_uni = X[:, idx_uni]
        x̄_uni = x_uni .- X̄
        ȳ_uni = y_uni .- Ȳ
        beta1_uni = (x̄_uni * x̄_uni') \ (x̄_uni * ȳ_uni)
        beta0_uni = Ȳ .- X̄'beta1_uni

        # opt
        # Random.seed!(1)
        e = Y .- Z'beta_uni
        # e = Y .- Z'beta_f
        # opt = sqrt.(sum(Z.^2, dims=1)) .* abs.(e')
        # opt = sqrt.(sum(((Z*Z') \ Z).^2, dims=1)) .* abs.(e')
        opt = abs.(e') .* sqrt.(lv) 
        # opt = opt ./ sum(opt)
        W_opt = 1 ./ opt; W_opt = W_opt ./ sum(W_opt)

        idx_opt = wsample(1:n, vec(opt), r, replace=true)
        z_opt = Z[:, idx_opt]
        y_opt = Y[idx_opt]
        opt_opt = opt[:,idx_opt]
        beta_opt = ((z_opt./opt_opt) * z_opt') \ (z_opt./opt_opt) * y_opt 
        # beta_opt[1] = (Ȳ .- X̄'beta_opt[2:end])[1]

        x_opt = X[:, idx_opt]
        X̄_opt = sum(X .* W_opt, dims=2)
        Ȳ_opt = sum(Y .* W_opt')
        x̄_opt = x_opt .- X̄_opt
        ȳ_opt = y_opt .- Ȳ_opt
        beta1_opt = ((x̄_opt./opt_opt) * x̄_opt') \ (x̄_opt./opt_opt) * ȳ_opt
        beta0_opt = Ȳ_opt .- X̄_opt'beta1_opt
        # beta0_opt = Ȳ .- X̄'beta1_opt

        # beta_opt[1] = (Ȳ_opt .- X̄_opt'beta_opt[2:end])[1]

        # lev
        # Random.seed!(1)
        W_lv = 1 ./ lv; W_lv = W_lv ./ sum(W_lv)
        idx_lv = wsample(1:n, vec(lv), r, replace=true)
        z_lv = Z[:, idx_lv]
        y_lv = Y[idx_lv]
        lv_lv = lv[:,idx_lv]
        beta_lv = ((z_lv ./ lv_lv) * z_lv') \ (z_lv ./ lv_lv) * y_lv 
        # beta_lv[1] = (Ȳ .- X̄'beta_lv[2:end])[1] 

        x_lv = X[:, idx_lv]
        X̄_lv = sum(X .* W_lv, dims=2)
        Ȳ_lv = sum(Y .* W_lv')
        x̄_lv = x_lv .- X̄_lv
        ȳ_lv = y_lv .- Ȳ_lv
        beta1_lv = ((x̄_lv ./ lv_lv) * x̄_lv') \ (x̄_lv ./ lv_lv) * ȳ_lv
        beta0_lv = Ȳ_lv .- X̄_lv'beta1_lv
        # beta0_lv = Ȳ .- X̄'beta1_lv

        # beta_lv[1] = (Ȳ .- X̄'beta_lv[2:end])[1] 
        # beta_lv[1] = (Ȳ_lv .- X̄_lv'beta_lv[2:end])[1] 
        # beta_lv[1] = (Ȳ .- X̄'beta1_lv)[1] 


        # iboss
        idx_iboss = getidx(X', r)
        z_iboss = Z[:,idx_iboss]
        y_iboss = Y[idx_iboss]
        beta_iboss = (z_iboss * z_iboss') \ (z_iboss * y_iboss)
        beta_iboss[1] = (Ȳ .- X̄'beta_iboss[2:end])[1]

        x_iboss = X[:, idx_iboss]
        x̄_iboss = x_iboss .- X̄
        ȳ_iboss = y_iboss .- Ȳ
        beta1_iboss = (x̄_iboss * x̄_iboss') \ (x̄_iboss * ȳ_iboss)
        beta0_iboss = Ȳ .- X̄'beta1_iboss

        Betas[:,rr,1] = beta_uni
        Betas[:,rr,2] = [beta0_uni; beta1_uni]
        Betas[:,rr,3] = beta_iboss
        Betas[:,rr,4] = [beta0_iboss; beta1_iboss]
        Betas[:,rr,5] = beta_lv
        Betas[:,rr,6] = [beta0_lv; beta1_lv]
        Betas[:,rr,7] = beta_opt
        Betas[:,rr,8] = [beta0_opt; beta1_opt]
        Betas[:,rr,9] = beta_f
        # next!(pgs)
    end
    return Betas
end
