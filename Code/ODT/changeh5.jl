using HDF5,DelimitedFiles


s_modes = ["d+","d+0","r+","d-o","dfr","d-","rfr","u","r-ab"]
s_parameter = ["A","ω","Γ","ϕ"]
param_names = [ "$(p)_$(m)" for m in s_modes, p in s_parameter ]
param_names = vcat(param_names...)

function change_h5(dir_molecule::String, dir_measurement::String)

    rd = readdir(joinpath(dir_molecule,dir_measurement))
    h5 = rd[occursin.(".h5",rd)][1]
    file =joinpath(dir_molecule,dir_measurement,h5)
    a = readdlm(joinpath(dir_molecule,dir_measurement,rd[occursin.("attenuation",rd)][1]))
    hold = readdlm(joinpath(dir_molecule,dir_measurement,rd[occursin.("hold",rd)][1]))[:] .|> Bool
    ref = readdlm(joinpath(dir_molecule,dir_measurement,rd[occursin.("ref",rd)][1]))[:]
    coha = try readdlm(joinpath(dir_molecule,dir_measurement,rd[occursin.("coha",rd)][1]))[:] catch end

    h5open(file,"r+") do io
        newgroup = create_group(io,"Fitting")
        
        newgroup["Attenuation"] = a
        newgroup["Parameter: Ref. Spectrum"] = ref
        attrs(newgroup)["Fitted Modes"] = s_modes
        attrs(newgroup)["Parameter Names: Ref. Spectrum"] = param_names
        attrs(newgroup)["Hold Modes"] = s_modes[hold]
        if coha !== nothing
            newgroup["Parameter: Coherent Artifact"] = coha
        end
    end
    
end

dir_odt = "../Data/Octadecanethiol(ODT)/"
dir_udt = "../Data/Undecanethiol(UDT)/" 
dir_ht  = "../Data/Hexanethiol(HT)/" 

dir_dmin = "Delay Scan/dminus pumped/"
dir_rmin = "Delay Scan/rminus pumped/"
dir_rfr = "Delay Scan/rfr pumped/"

dir_wl = "Wavenumber Scan/"


#change_h5(dir_odt,dir_dmin)
#change_h5(dir_odt,dir_rmin)
#change_h5(dir_odt,dir_rfr)
#change_h5(dir_odt,dir_wl)
#change_h5(dir_udt,dir_dmin)
#change_h5(dir_udt,dir_rmin)
#change_h5(dir_udt,dir_rfr)
#change_h5(dir_ht,dir_dmin)
#change_h5(dir_ht,dir_rmin)
#change_h5(dir_ht,dir_rfr)


h5open(joinpath(dir_odt,dir_dmin,"DL-ODT-012-001.h5"), "r") do fid read(fid) end