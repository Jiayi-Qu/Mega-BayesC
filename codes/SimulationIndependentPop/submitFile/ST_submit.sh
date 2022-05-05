currentpath="/home/MegaLMM_Bayesian/BayesC_Simulation/"
Jcode="$currentpath/Codes/ST_analysis.jl"
nFactors=(2 6 9)
nTraitsPerFactor=(2 20)
n_qtl=(10 20 30)

for i in "${nFactors[@]}"
do
	echo "nFactors = $i"
	for t in "${nTraitsPerFactor[@]}" 
	do
		nTraits=$((i * t))
		echo "nTraits = $nTraits"
		for j in "${n_qtl[@]}"
		do
			echo "snp number = $j"
			data_path="$currentpath/${i}factor_${nTraits}trait_${j}snp_0.95h2F"
			for r in {1..10}
			do
				echo "repeat $r"
				echo "#SBATCH -J snp${j}_rep$r" > job_name
				echo "#SBATCH -D $data_path" > path_name
				echo module load julia/1.3.0 >load_module
				echo julia $Jcode $data_path $nTraits $i $j $r  > ender
				file="ST_trait${nTraits}_factor${i}_snp${j}_rep$r.sh"
				cat header job_name path_name load_module ender > $file
				chmod u+x $file
				sbatch $file
			done
		done
	done
done
