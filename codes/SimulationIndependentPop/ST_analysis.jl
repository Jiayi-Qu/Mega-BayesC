#using Pkg
#Pkg.add("StatsBase")
using JWAS,DelimitedFiles,CSV,DataFrames,Random,StatsBase,Statistics;

data_path =  ARGS[1]
nTraits = ARGS[2]
nFactors = ARGS[3]
numberofsnp = ARGS[4]
replicate = ARGS[5]

println("nFactor = $nFactors")
println("nTrait = $nTraits")
println("nSNP = $numberofsnp")
println("replicate = $replicate")

PHENO_ORIGIN = CSV.read(data_path*"/Y_"*nFactors*"Factors_"*nTraits*"Traits_rep$replicate.csv",delim = ',',header=true);
column_names = names(PHENO_ORIGIN)
rename!(PHENO_ORIGIN, :Column1 => :ID);

try
    mkdir(data_path*"/ST_rep$replicate")
catch
    "a folder has been existed."
end

cd(data_path*"/ST_rep$replicate")

model_equation = "focal = intercept"
model = build_model(model_equation)
genofile = data_path*"/X2_rep$replicate.csv"
add_genotypes(model,genofile,separator=',',header=true)
outputMCMCsamples(model)
out = runMCMC(model, PHENO_ORIGIN, Pi=0.95,estimatePi=true,methods="BayesC",
    chain_length=10000,outputEBV=true,output_samples_file=data_path*"/ST_rep$replicate/"); 
