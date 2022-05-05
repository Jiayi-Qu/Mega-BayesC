using DataFrames,CSV,JWAS,Statistics,DelimitedFiles,Pkg,Random,Distributions,LinearAlgebra;

analysis_path = ARGS[1]
data_path = ARGS[2]


geno_file = data_path * "geno_file_LDclumped.csv"

pheno = CSV.read(data_path * "BayesMega_FT10.csv")
rename!(pheno, :phenotype_value => :FT)

pheno = pheno[:, [:accession_id, :FT]]

try
mkdir(analysis_path*"STBayesC")
catch
    "a folder has been existed."
end

Random.seed!(123);
cd(analysis_path*"STBayesC")

# run BayesC
model = build_model("FT = intercept")
add_genotypes(model,geno_file,header=true,separator=',');

@time output=runMCMC(outputEBV=true,model,pheno, 
    chain_length=50_000,burnin=10_000, estimatePi = true, #Pi = 0.99, output_samples_file=analysis_path*"replicate$rep/Day$i/",
    methods="BayesC",double_precision=true, 
    output_samples_file=analysis_path*"STBayesC/");
