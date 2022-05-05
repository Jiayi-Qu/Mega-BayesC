using DataFrames,CSV,JWAS,Statistics,DelimitedFiles,Pkg,Random,Distributions,LinearAlgebra;

analysis_path = ARGS[1]
data_path = ARGS[2]
pop = ARGS[3]
pop = parse(Int, pop)

## LDclumped 
println("LDclumped")

geno_file = data_path * "replicate$pop/geno_file_LDclumped_rep$pop.csv"

pheno = CSV.read(data_path * "replicate$pop/FT10_sim_gexp_ind_rep$pop.csv")
rename!(pheno, :phenotype_value => :FT)

pheno = pheno[:, [:accession_id, :FT]]

try
mkdir(analysis_path*"LDclumped/replicate$pop/STBayesC")
catch
    "a folder has been existed."
end

Random.seed!(123);
cd(analysis_path*"LDclumped/replicate$pop/STBayesC")

# run BayesC
model = build_model("FT = intercept")
add_genotypes(model,geno_file,header=true,separator=',');

@time output=runMCMC(outputEBV=true,model,pheno, 
    chain_length=50_000,burnin=10_000, estimatePi = true, 
    methods="BayesC",double_precision=true, 
    output_samples_file=analysis_path*"LDclumped/replicate$pop/STBayesC/");


