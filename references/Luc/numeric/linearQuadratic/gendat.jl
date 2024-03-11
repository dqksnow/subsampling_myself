function gendat(n::Int, case::Int=1; beta0::Vector=ones(10))
    d = length(beta0)
    Zt = rand(n)
    Z  = [Zt Zt.^2]
    Y = beta0[1] .+ Z * beta0[2:end] + randn(n)
    return Z, Y
end
