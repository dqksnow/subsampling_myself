using DelimitedFiles

for d in [2, 20]
    tb = Matrix{Any}(undef, 7, 0)
    for i in 1:3
        rs = readdlm("output/d$(d)case$(i).csv")[1:end-3,:]
        tb = [tb rs]
    end
    fl = typeof.(tb) .== Float64
    tb[fl] .*= 1000
    writedlm("output/00d$(d)res.csv", tb)
end


for d in [2, 20]
    tb = Matrix{Any}(undef, 6, 0)
    nams = readdlm("output/d$(d)case1.csv")[2:end-3,1]
    for i in 1:3
        rs = readdlm("output/d$(d)case$(i).csv")
        tb = [tb rs[2:end-3,2:end]]
    end
    writedlm("output/00Td$(d)res.csv", [permutedims(nams); 100tb'])
end
