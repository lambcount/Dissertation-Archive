# Dissertation-Archive

This repository contains the data and models presented in my dissertation (accessible via [DuEPublico link]).
The entire repository is also archived on CERN's cloud storage (DOI: ).

All relevant parameters and spectra data can be found in HDF5 files under the `Data/` directory, alongside the corresponding fitting results discussed in the dissertation. Parameters relevant for reproducing these fits are also included in these files.

The source code provided under the `Code/` directory was developed using Julia version 1.9.4.

Below you can find a brief tutorial detailing how to reproduce the presented analyses and models.

## Tutorial

### Activate the Julia Environment

All Julia packages used for this project are stored in the provided `Manifest.toml.`
Activate the Julia environment by navigating into the `Code/` folder and running:

```julia
    julia> using Pkg
    julia> Pkg.activate(".")
```

If you initialize this environment for the first time, you must install all required packages using:

```julia
    julia> using Pkg
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()
```

### Loading Experiment Data

All experimental data is stored in the `Data/` directory. The data is provided as HDF5 files, which can be loaded as follows:

```julia
    julia> using HDF5
    julia> h5open(dir, "r") do fid read(fid) end
    Dict{String, Any} with 2 entries:
        "Fitting" => Dict{String, Any}
        "Data"    => Dict{String, Any}
```

Each HDF5 file contains the spectral data (`"Data"`) and the corresponding fitting results (`"Fitting"`). The fitting results were obtained using the `fit_signal_const_ref_attenuation` function provided in `Code/SFG.jl`.

In some fits, certain attenuation parameters were held constant (not varied); these are documented under attributes of `"Fitting"`.

Example workflow to load the data and fit the signal with the `fit_signal_const_ref_attenuation` function:

```julia
    julia> using HDF5,LsqFit
    julia> data = h5open(dir, "r") do fid read(fid)["Data"] end 
    julia> ref_spectrum = h5open(dir, "r") do fid read(fid)["Fitting"]["Parameter: Ref. Spectrum"] end 
    julia> hold_modes =  h5open(dir, "r") do fid HDF5.attributes(fid["Data"])["Hold Modes"] |> read end
    julia> modes = h5open(dir, "r") do fid HDF5.attributes(fid["Data"])["Fitted Modes"] |> read end
    julia> hold = [x in hold_modes ? 1 : 0 for x in modes]
    julia> fitted_data = fit_signal_const_ref_attenuation(ref_spectrum,data,hold=hold)
```

The attenuations extracted from `fitted_data` are already stored under `["Fitting"]["Attenuation"]` within the respective HDF5 files.

### Models

The Julia code used for the models described in my dissertation is provided in `Code/Models.jl`. These models were optimized using the global optimization algorithm `:adaptive_de_rand_1_bin_radiuslimited` provided by the `BlackBoxOptim.jl` package.

Below is an example demonstrating the fitting of `Model 3` to previously determined attenuation data. In the parameter `sel_$(experiment)`, modes with continuous behavior and no signal enhancement (as discussed in the dissertation) are selected. The indices refer to the specific modes used in the model function.
> ⚠️  **Important:** Since the modes `dplus` and `dplusomega` do not exhibit continuous behavior, the index `[1]` refers to the `rplus` mode, not the `dplus` mode!

The parameter `coha` is an array containing two elements: `coha[1]` is the index of the coherent mode, and `coha[2]` indicates its respective time indices. This parameter was only used for the rfr-pumped experiment, relevant specifically for the molecule ODT.

Example usage:

```julia
    julia> using BlackBoxOptim, HDF5
    julia> include("Models.jl")
    julia> exp_dmin = h5open(dir_ODT_dmin, "r") do fid read(fid) end
    julia> exp_rfr = h5open(dir_ODT_rfr, "r") do fid read(fid) end
    julia> exp_rmin = h5open(dir_ODT_rmin, "r") do fid read(fid) end
    julia> coha = [3,[19:31]]       # Mode [3] is the coherent mode, and [19:31] are the time indices
    julia> dt = 0.1                 # in ps
    julia> sel_dmin = [1,5,6,7]     # selected modes for the dmin pumped experiment
    julia> sel_rfr = [1,3,5,7]      # selected modes for the rfr pumped experiment
    julia> sel_rmin = [1,2,5,6,7]   # selected modes for the rmin pumped experiment

    julia> model3_3_dmin,model3_3_rfr,model3_3_rmin = blackbox_model_3_3(
        exp_dmin["Data"]["dltime"],exp_dmin["Fitting"]["Attenuation"][sel_dmin,:],sel_dmin,
        exp_rfr["Data"]["dltime"],exp_rfr["Fitting"]["Attenuation"][sel_rfr,:],sel_rfr,
        exp_rmin["Data"]["dltime"],exp_rmin["Fitting"]["Attenuation"][sel_rmin,:],sel_rmin,
        coha,
        dt=dt
    )
```

The fitted models are in the `Data/$(molecule)/Model` directory, where `molecule` is the name of the molecule used in the experiment. The models are stored in JLD2 File to preserve the data structure in `Code/Reservoir Model Struct.jl`. There is a minimal working example on how to load the models in `Code/Load Model Example.jl`.

In addition to the provided fitting examples, it is also possible to define custom models using the general model functions included in `Models.jl`. You can adjust these model functions with user-defined parameters to suit your individual requirements. Detailed information on available parameters and usage instructions can be found directly within the comments and documentation of each function inside the `Models.jl` file. With the kwargument `mode`, you can choose which mode to pump. The default value is `mode = :dmin`. The kwargument `dt` allows you to set the time step for the model. The default value is `dt = 0.1` ps.

```julia
    julia> include("Models.jl")
    julia> time = collect(-30:0.1:300)  # time vector in ps
    julia> parameters = rand(43)        # random parameters for the model
    julia> model = model3_3(            # construct a model where rplus mode is pumped
        time, 
        parameters, 
        mode=:rplus, 
        dt = 0.1
        ) 
```

> ⚠️  **Important:** Please replace all placeholders (such as `dir_ODT_dmin`, `dir_ODT_rfr`, and `dir_ODT_rmin`) with the actual paths and filenames in `Data`.
