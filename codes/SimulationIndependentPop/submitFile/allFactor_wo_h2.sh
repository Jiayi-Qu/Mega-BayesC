currentpath="/home/MegaLMM_Bayesian/BayesC_Simulation/"
Rcode="$currentpath/Codes/simulation_allFactor_wo_h2.R"
code_path="$currentpath/Codes"
nFactors=(2 6 9)
nTraitsPerFactor=(2 20)
n_qtl=(10 20 30)
#snp_per_window=20
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
			simulation_path="$currentpath/allFactor_wo_h2/${i}factor_${nTraits}trait_${j}snp_0.95h2F"
			mkdir $simulation_path
			for r in {1..10}
			do
				echo "repeat $r"
				echo "#SBATCH -J snp${j}_rep$r" > job_name
				echo "#SBATCH -D $currentpath/allFactor_wo_h2/" > path_name
				echo module load R >load_module
				echo Rscript $Rcode $simulation_path $code_path $nTraits $i $j $r  > ender
				file="allFactor_woh2_trait${nTraits}_factor${i}_snp${j}_rep$r.sh"
				cat header job_name path_name load_module ender > $file
				chmod u+x $file
				sbatch $file
			done
		done
	done
done
