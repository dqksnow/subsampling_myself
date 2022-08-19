using Plots, DelimitedFiles, LaTeXStrings
# Plots.scalefontsizes(1.8)
# cd("..")

# ncase = [1:5; 10][1:5]
# n=[fill(10^5, 4); 4188261][1:5]
# # Ms = [1, 2, 5, 10, 100]
# # label = ["dc: M=" .* map(string, Ms); "uniform"]
nmd = 5
N = 10^5
ncase = 1:nmd
Nss=fill(N, nmd)
label = [L"\mathbf{Sequential}" L"\mathbf{Uniform}" L"\mathbf{IBOSS}" L"\mathbf{OSMAC}" L"\mathbf{OSMAC-Rep}"]
for (i, case) in enumerate(ncase)
    rs = readdlm("output/case$(case).csv")
    pl = plot(rs[:,1]/Nss[i], log.(rs[:,2:end-1]),
              label=label, lw=3, m=(10,:auto),
              tickfontsize=18, xguidefontsize=20,
              legendfontsize=15, grid=false# ,
              # legend= case ==10 ? true : false
              )
    xlabel!("\\alpha")
    ylabel!("log(MSE)")
    savefig(pl, "output/0case$(case).pdf")
end

# rec = fill(NaN, nmd, length(nss))
# # @time cal = simu(X, Y, nss, rpt, nmd)
# for m in 1:nmd, (idn, n) in enumerate(nss)
#     fname = "output/csv/case$(case)N$(N)method$(m)n$(n).csv"
#     betas = readdlm(fname)
#     rec[m,idn] = sum(mean((betas .- beta0).^2, dims=2))
#     # readdlm(fname, Betas[:,idn,:,m])
# end
# [nss rec']

# lg = plot(ones(1, ncase+1), label=label, lw=3, m=(10,:auto),
#           showaxis=false, grid=false,
#           legend=:inside, legendfontsize=18)
# scatter!(ones(1, 1), color=:white,
#          markerstrokecolor=:white, m=(20), label="")
# savefig(lg, "output/case$(ncase+1).pdf")



pdfs = "output/0case" .* map(string, ncase) .* ".pdf"
run(`pdftk $pdfs output output/0mse_linear.pdf`)
# run(`pdftk $pdfs output $(homedir())/Dropbox/work/figures/mse_linear.pdf`)

#
# rs = readdlm("output/case1.csv")
