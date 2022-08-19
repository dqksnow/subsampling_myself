# using LinearAlgebra
# using Distributions
# using Random
#
# n = 10^6
# beta0 = beta_0  = ones(100)
# d = p = size(beta_0)[1]
# corr  = 0.5
# sigmax = [corr+(1-corr)*(i==j) for i in 1:d, j in 1:d]
# Random.seed!(1);
#
# k = 1000
# # genX = MvNormal(zeros(d), sigmax);
# Z = rand(MvNormal(zeros(d), sigmax), n);
# Z = convert(Array{Float64,2}, Z');
# # Z = X;
# X = [ones(n) Z];
# beta0 = [1; beta_0];
# p = d + 1
# Y  = vec(X * beta0) + rand(Normal(0, 3), n);
#
# # Full
# @time sum(inv(X' * X) * (X' * Y));

function getidxo(x, k)
    (n,p) = size(x)
    idx = Int[]
    r = cld(k, 2p)
    tmp2 = x[:,1]
    l = sort!(tmp2; alg=PartialQuickSort(r))[r]
    u = sort!(tmp2; alg=PartialQuickSort(r), rev=true)[r]
    bl = 1; bu = 1;
    for s in 1:n
        if bl > r && bu > r
            break
        elseif bl <= r && x[s,1] <= l
            push!(idx, s)
            bl += 1
        elseif bu <= r && x[s,1] >= u
            push!(idx, s)
            bu += 1
        end
    end
    for j in 2:p
        sort!(idx)
        tl = n - length(idx)
        tmp = Vector{Float64}(undef, tl)
        t = 1; v = 1;
        for s in 1:n
            # global v, t
            if s != idx[v]
                tmp[t] = x[s,j]
                t += 1
            elseif v < length(idx)
                v += 1
            end
        end
        l = sort!(tmp; alg=PartialQuickSort(r))[r]
        u = sort!(tmp; alg=PartialQuickSort(r), rev=true)[r]
        v = 1; bl = 1; bu = 1;
        for s in 1:n
            if length(idx) >= k || (bl > r && bu > r)
                break
            elseif v <= length(idx) && s == idx[v]
                v += 1
            elseif bl <= r && x[s,j] <= l
                push!(idx, s)
                bl += 1
            elseif bu <= r && x[s,j] >= u
                push!(idx, s)
                bu += 1
            end
        end
    end
    return idx
end
# @time old = sum(getidxo(Z, 10000))
# x = Z[:,2:end]
# k = k2
# getidxo(Z[:,2:end], k2)

function getidx(x, k)
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

#
# function getidx(x::Matrix, k::Int64)
#     (n,p) = size(x)
#     # d = p - 1
#     idx = Array{Int64}(undef, k)
#     counter = 1
#     r = cld(k, 2(p-1))
#     tmp2 = x[:,2]
#     l = sort!(tmp2; alg=PartialQuickSort(r))[r]
#     u = sort!(tmp2; alg=PartialQuickSort(r), rev=true)[r]
#     bl = 1; bu = 1;
#     for s in 1:n
#         # global counter
#         if bl > r && bu > r
#             break
#         elseif bl <= r && x[s,2] <= l
#             idx[counter] = s
#             counter += 1
#             bl += 1
#         elseif bu <= r && x[s,2] >= u
#             idx[counter] = s
#             counter += 1
#             bu += 1
#         end
#     end
#     for j in 3:p
#         tl = n - counter + 1
#         tmp = Vector{Float64}(undef, tl)
#         t = 1; v = 1;
#         for s in 1:n
#             # global v, t
#             if s != idx[v]
#                 tmp[t] = x[s,j]
#                 t += 1
#             elseif v < counter - 1
#                 v += 1
#             end
#         end
#         l = sort!(tmp; alg=PartialQuickSort(r))[r]
#         u = sort!(tmp; alg=PartialQuickSort(r), rev=true)[r]
#         v = 1; bl = 1; bu = 1;
#         for s in 1:n
#             # # global v, k, counter
#             # if counter > k
#             #     break
#             # elseif s == idx[v]
#             #     v += 1
#             # elseif x[s,j] <= l || x[s,j] .>= u
#             #     idx[v+1:counter] = idx[v:counter-1]
#             #     idx[v] = s
#             #     v += 1
#             #     counter += 1
#             # end
#             if counter > k || (bl > r && bu > r)
#                 break
#             elseif s == idx[v]
#                 v += 1
#             elseif bl <= r && x[s,j] <= l
#                 idx[v+1:counter] = idx[v:counter-1]
#                 idx[v] = s
#                 v += 1
#                 counter += 1
#                 bl += 1
#             elseif bu <= r && x[s,j] >= u
#                 idx[v+1:counter] = idx[v:counter-1]
#                 idx[v] = s
#                 v += 1
#                 counter += 1
#                 bu += 1
#             end
#         end
#     end
#     return idx
# end
# @time new = sum(getidx(X, 1000));
#
#
# # @time tmp = [X[i,2] for i in 1:n];
# #
# # tmp = []
# # @time for i in 1:n
# #     push!(tmp, X[i,2])
# # end;
# #
# # @time tmp = Vector{Float64}(undef,n);
# # @time for i in 1:n
# #     tmp[i] = X[i,2]
# # end;
# #
# #
# # ccall((:clock, ), Int32, ())
#
# # function oT(x::Matrix, y::Vector)
# #     k = 1000; rpt = 100;
# #     (n,p) = size(x)
# #     norm = Array(Float64, n);
# #     for rr in 1:rpt
# #         for i in 1:n
# #             norm[i] = 0
# #             for j in 1:p
# #                 norm[i] = norm[i] + x[i,j]^2
# #             end
# #             norm[i] = sqrt(norm[i])
# #         end
# #         mk = select!(norm, k, rev=true)
# #         idx_oT = find(norm .>= mk)
# #         x_oT = x[idx_oT,:]
# #         y_oT = x[idx_oT]
# #         estimate = inv(x_oT' * x_oT) * (x_oT' * y_oT)
# #     end
# #     return estimate
# # end
# # @time oT(X, Y)
# #
# #
# #
# #
# # function oD(x::Matrix, y::Vector)
# #     k = 1000; d = 10; p = d +1; rpt = 100;
# #     for rr in 1:rpt
# #     z = copy(x)
# #         idx = 0
# #         r = 10 ## convert(Int, k / d / 2)
# #         mru = select!(z[:,2], r, rev=true)
# #         mrl = select!(z[:,2], r)
# #         idx = find(x[:,2] .>= mru)
# #         idx = [idx, find(x[:,2] .<= mrl)]
# #         for j in 3:p
# #             mru = select!(z[:,j], r, rev=true)
# #             mrl = select!(z[:,j], r)
# #             idx = [idx, find(x[:,j] .>= mru)]
# #             idx = [idx, find(x[:,j] .<= mrl)]
# #         end
# #         ## idx = unique(idx)
# #         ## x_oD = x[idx,:]
# #         ## y_oD = y[idx]
# #         ## estimate = inv(x_oD' * x_oD) * (x_oD' * y_oD)
# #     end
# #     # return 1
# # end
# # @time oD(X, Y)
# #
# #
# #
# # x = rand(1:500, 10)
# # k = 5
# # k2 = 5:10
# # s = sort(x; alg=QuickSort)
# # ps = sort(x; alg=PartialQuickSort(k))
# # qs = sort(x; alg=PartialQuickSort(k2))
#
# # m = select(x, k, rev=true)
# # find(x .>= m)
