using Pkg 
# Activate the project environment
# Make sure to run this script, when pwd() is either "Code" or "Disseration-Archive"
if Sys.iswindows() 
    _pwd = split(pwd(),"\\")[end]
    active_dir = split(Base.active_project(),"\\")[end-1]
    if _pwd != "Code"
        cd("Code")
        _pwd = split(pwd(),"\\")[end]
        if _pwd != "Code"
            error("Wrong directory")
        else
            active_dir = split(Base.active_project(),"\\")[end-1]
        end
    end
    if active_dir != "Code"
        Pkg.activate(".")
    end
    elseif Sys.isapple()
        _pwd = split(pwd(),"/")[end]
        active_dir = split(Base.active_project(),"/")[end-1]
        if _pwd != "Code"
            cd("Code")
            _pwd = split(pwd(),"/")[end]
            if _pwd != "Code"
                error("Wrong directory")
            else
                active_dir = split(Base.active_project(),"/")[end-1]
            end
        end
        if active_dir != "Code"
            Pkg.activate(".")
        end
    end

using JLD2,FileIO,DataFrames
include("Reservoir Model Struct.jl")

# Load the model from the JLD2 file
molecule = "Octadecanethiol(ODT)"
dir_molecule = "../Data/$molecule/"

dir_model = joinpath(dir_molecule,"Model","models.jld2")
model = load(dir_model)["models_dir"]

# Extract the different models variants
k_variants = keys(model) |> collect
k_variant = k_variants[1]  

#Extract Pump-Probe-Delay, Bleach, Population and used parameters
# Bleach and Population are mxn matrices, where m corrensponds to the model pump-probe delay and n to different modes in the model. 
#See model function for more details.
bleach  = model[k_variant]["$k_variant dmin pumped"]
population = model[k_variant]["$k_variant dmin pumped"].Population
pump_probe_delay = model[k_variant]["$k_variant dmin pumped"].Time
parameters = model[k_variant]["$k_variant dmin pumped"].Parameter

# Extract Dictionary of the model
model_dict = model[k_variant]["$k_variant dmin pumped"].Dictionary