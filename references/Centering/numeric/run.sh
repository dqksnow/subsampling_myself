#! /bin/bash

[ ! -d "./output" ] && mkdir output

for d in 2 20
do
		echo "dims: $d"
		export dims=$d
		for i in {1..3} # 10
		do
				echo "case: $i"
				export case=$i
				nohup julia -t 15 run.jl > output/res${d}-${i}.out &
				sleep 0.1s
		done
done

# echo "output 00mse.pdf" >> output/00combinePDF.sh
# echo "sleep 0.1s" >> output/00combinePDF.sh
# # echo "cp 00mse.pdf ~/Dropbox/work/Ma/revision1/figures/00mse.pdf" >> output/00combinePDF.sh
# chmod +x output/00combinePDF.sh

