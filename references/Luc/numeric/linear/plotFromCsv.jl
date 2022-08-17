initime = time()
using DelimitedFiles, Statistics, LinearAlgebra
using Plots, LaTeXStrings
# Threads.nthreads()
BLAS.set_num_threads(4)
# ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())

# Plots.scalefontsizes(1.8)
# cd("..")

# ncase = [1:5; 10][1:5]
# Nss=[fill(10^5, 4); 4188261][1:5]
# # Ms = [1, 2, 5, 10, 100]
# # label = ["dc: M=" .* map(string, Ms); "uniform"]
N = 10^5
d = 10
beta0 = ones(d)
nmd = 5
ncase = 1:6
nss = [2000, 5000, 10000, 20000, 50000]
beta0 = ones(d)

Nss=fill(N, 6)
label = [L"\mathbf{Sequential}" L"\mathbf{Uniform}" L"\mathbf{IBOSS}" #=
         =# L"\mathbf{OSMAC}" L"\mathbf{OSMAC-Rep}"]
cmse = ["slope", "intercept", "mse"][3]

# Threads.@threads for i in 1:length(ncase)
for i in 1:length(ncase)
    case = ncase[i]
    rec = fill(NaN, nmd, length(nss))
# @time cal = simu(X, Y, nss, rpt, nmd)
    for m in 1:nmd, (idn, n) in enumerate(nss)
        fname = "output/csv/case$(case)N$(N)method$(m)n$(n).csv"
        betas = readdlm(fname)
        if cmse == "slope"
            rec[m,idn] = sum(mean((betas .- beta0).^2, dims=2)[2:end])
        elseif cmse == "intercept"
            rec[m,idn] = sum(mean((betas .- beta0).^2, dims=2)[1])
        elseif cmse == "mse"
            rec[m,idn] = sum(mean((betas .- beta0).^2, dims=2))
        end
    end
    rs = [nss rec']
    # rs = readdlm("output/case$(case).csv")
    pl = plot(rs[:,1]/Nss[i], log.(rs[:,2:end-2]),
              label=label, lw=3, m=(10,:auto),
              tickfontsize=18, xguidefontsize=20,
              legendfontsize=15, grid=false# ,
              # legend= case ==10 ? true : false
              )
    xlabel!("\\alpha")
    ylabel!("log(MSE)")
    savefig(pl, "output/0case$(case).pdf")
end


# lg = plot(ones(1, ncase+1), label=label, lw=3, m=(10,:auto),
#           showaxis=false, grid=false,
#           legend=:inside, legendfontsize=18)
# scatter!(ones(1, 1), color=:white,
#          markerstrokecolor=:white, m=(20), label="")
# savefig(lg, "output/case$(ncase+1).pdf")

pdfs = "output/0case" .* map(string, ncase) .* ".pdf"
run(`pdftk $pdfs output output/0mse_linear_$(cmse).pdf`)
run(`cp output/0mse_linear_$(cmse).pdf ../draft/figures/`)


# run(`pdftk $pdfs output $(homedir())/Dropbox/work/figures/mse_linear.pdf`)

#
# rs = readdlm("output/case1.csv")

# initime = time()
# # Threads.@threads 
# for i in 1:length(ncase)
#     case = ncase[i]
#     rec = fill(NaN, nmd, length(nss))
#     for m in 1:nmd, (idn, n) in enumerate(nss)
#         fname = "output/csv/case$(case)N$(N)method$(m)n$(n).csv"
#         betas = readdlm(fname)
#         rec[m,idn] = sum(mean((betas .- beta0).^2, dims=2)[2:end])
#     end
#     rs = [nss rec']
#     # rs = readdlm("output/case$(case).csv")
#     pl = plot(rs[:,1]/n[i], log.(rs[:,2:end-1]),
#               label=label, lw=3, m=(10,:auto),
#               tickfontsize=18, xguidefontsize=20,
#               legendfontsize=15, grid=false# ,
#               # legend= case ==10 ? true : false
#               )
#     xlabel!("\\alpha")
#     ylabel!("log(MSE)")
#     savefig(pl, "output/0case$(case).pdf")
# end
print("\n Total time used: $(round(time() - initime)) seconds")
