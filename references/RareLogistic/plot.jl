using Plots, DelimitedFiles, LaTeXStrings
# Plots.scalefontsizes(1.8)

N = 10^6 ÷ 2
ncase = 1:4
# # label = [L"\mathbf{Sequential}", L"\mathbf{Uniform}",
# # L"\mathbf{OSMAC}", L"\mathbf{OSMAC-Rep}"]
label = ["uniW" "uniLik" "optW" "optLik" "LCC" "Full"]
tp = 1:5

for (i, case) in enumerate(ncase)
    rs = readdlm("output/case$(case).csv")# [1:end-1,:]
    rs[:,2:end] = log.(rs[:,2:end])
    plc = plot(rs[:,1]./N, rs[:,tp.+1], # size = 1.2 .*(400, 600),
              label=label[:,tp], lw=3, m=(9,:auto),
              tickfontsize=16, xguidefontsize=18, yguidefontsize=18,
              legendfontsize=14, grid=false, thickness_scaling=1,
              # xlabel="sampling rate",
               ylabel="log(MSE)"# , legend=:topright
               ,legend= case == 1 ? :topright : false# ,
               # title = "Correct Pilot", top_margin=1
               )
    annotate!(0.01, sum(extrema(rs[:,tp.+1]) .* [0.05, 0.95]),
              text("Consistent Pilot", :center, 16))
    fullmse = readdlm("output/full-case$(case).csv")
    hline!(log(fullmse), label="Full")
    savefig(plc, "output/0case$(case).pdf")
end

# lg = plot(ones(1, ncase+1), label=label, lw=3, m=(10,:auto),
#           showaxis=false, grid=false,
#           legend=:inside, legendfontsize=18)
# scatter!(ones(1, 1), color=:white,
#          markerstrokecolor=:white, m=(20), label="")
# savefig(lg, "output/case$(ncase+1).pdf")

pdfs = "output/0case" .* map(string, ncase) .* ".pdf"
onefile = "00mse.pdf"
run(`pdftk $pdfs output output/$(onefile)`)
paperdir = "$(homedir())/Dropbox/Apps/Overleaf/MSLE-Binary/figures"
run(`cp output/$(onefile) $(paperdir)/$(onefile)`)

