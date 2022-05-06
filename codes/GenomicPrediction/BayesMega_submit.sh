code="MegaLMM_Krause_geno_BayesC.R"
for i in {1..20}
do
	echo "replicate $i"
        echo "#SBATCH -J MegaBC$i" > job_name
	echo "#SBATCH -D /home/carol666/MegaLMM_Bayesian/Krause/Horseshoe/" > path_name
	echo module load R >load_module
        echo Rscript $code $i > ender
        file="MegaBayesCn$i.sh"
        cat header job_name path_name load_module ender > $file
       chmod u+x $file
       sbatch $file
done


