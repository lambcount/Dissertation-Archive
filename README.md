# Dissertation-Archive

    Im Folgenden werden die Daten und Modelle, welche ich, meiner Arbeit von mir (DuEPublico Link) vorgestellt habe, archiviert. 
    Das Golgende Repository ist auch auch der Cern eigenen Cloud bereitgestellt. (DOI:)
    Alle relevanten Parameter meiner Spektren sind die HDF5 Dateien unter Data-> zu finden. Dort sind ebenfalls die in der Arbeit dargestellten Fits 
    bereitgestellt. Die relevanten parameter f"uer die Anfertigung dieser sind dort ebvenfalls enthalten.
    Der in Code -> enthaltene Programmcode wurde mit der Julia Version 1.9.4 angefertigt.
    Im Folgenden gebe ich eine kurzes Tutorial wie einzelnen Darstellungen angefertigt wurden.

## Tutorial

### Aktiviere das Julia-Eviroment

    Alle f"ur diese Arbeit verwendeten Pakete sind Maninfest.toml file abgespeichert. Das Enviroment
    kann folgenderma"sen aktiviert werden. In Julia in den Ordner "Code" navigieren. Dann:

```julia
    julia> using Pkg
    julia> Pkg.activate(".")
```

    Wird das Enviroment das 1. Mal initilasiert so m"ussen alle erforderlichen Pakete zun"achst 
    installiert werden. Dazu:

```julia
    julia> using Pkg
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()
```

### Experimente Laden

    Die durchgef"uhrten Experiment sind unter Data-> zu finden. Alle Experimente sind als HDF5 Dateien 
    gespeichert und k"onnen wie folgt geladen werden.

```julia
    julia> using HDF5
    julia> h5open(dir, "r") do fid read(fid) end
    Dict{String, Any} with 2 entries:
        "Fitting" => Dict{String, Any}
        "Data"    => Dict{String, Any}
```

   In dieser Datei sind sowohl die Spektren Informationen ("Data") als auch die Anpassungen ("Fitting").

   Die Anpassungen wurden mithilfe der Funktion "fit_signal_const_ref_attenuation" erhalten. Welche
   in Code->SFG.jl enthalten ist. F"ur manche Anpassungen wurden verschiedenene 
   Abschw"achungen nicht varriert. Diese sind in der Datei unter den Attributen von "Fitting"->" zu finden.

```julia
    julia> using HDF5,LsqFit
    julia> data = h5open(dir, "r") do fid read(fid)["Data"] end 
    julia> ref_spectrum = h5open(dir, "r") do fid read(fid)["Fitting"]["Parameter: Ref. Spectrum"] end 
    julia> hold_modes =  h5open(dir, "r") do fid HDF5.attributes(fid["Data"])["Hold Modes"] |> read end
    julia> modes = h5open(dir, "r") do fid HDF5.attributes(fid["Data"])["Fitted Modes"] |> read end
    julia> hold = [x in hold_modes ? 1 : 0 for x in modes]
    julia> fitted_data = fit_signal_const_ref_attenuation(ref_spectrum,data,hold=hold)
```

    Aus fitted_data k"onnen dann die Abschw"achungen extrahiert werden. Diese sind jedoch bereits unter Fitting->Attenuation gespeichert.

### Modelle 

    F"ur die in meiner Arbeit verwendeten Modelle wurde der entsprechende Code in Code->Models.jl verwendet. 
    Diese Modelle wurden mit einem globalen Optimierer aus dem BlackBoxOptim.jl Paket optimiert. 
    Im Folgenden wird ein Beispiel f"ur den Fit des Modell 3 an die zuvor bestimmten Absch"w"achungen gegeben. Unter "sel" wurden die in der Arbeit er"wahnten Transienten der Moden ausge"w"ahlt, welche keine Signalerh"ohung und einen kontinuierlichen Verlauf aufweisen.  Diese werden als Indice der verwendeten Moden angegeben. Wichtig!: Da die Moden d+ und d+0 keine kontinulieren Verl"aufe bildeten, steht der Index [1] f"ur die rplus Mode und nicht f"ur die d+ Mode!.
    Der Parameter "coha" ist ein Array mit 2 Eintr"agen. "coha[1]" ist der Index der Mode in und "coha[2]" die Indices in den Zeiten der jeweiligen Mode. Der wurde lediglich in der rfr-gepumpten Messung angewendet, und war nur f"ur das Molek"ul ODT relevant.

```julia
    julia> using BlackBoxOptim, HDF5
    julia> include("Models.jl")
    julia> exp_dmin = h5open(dir_ODT_dmin, "r") do fid read(fid) end
    julia> exp_rfr = h5open(dir_ODT_rfr, "r") do fid read(fid) end
    julia> exp_rmin = h5open(dir_ODT_rmin, "r") do fid read(fid) end
    julia> coha = [3,[19:31]]
    julia> dt = 0.1 # in ps
    julia> sel_dmin = [1,5,6,7]
    julia> sel_rfr = [1,3,5,7]
    julia> sel_rmin = [1,2,5,6,7]

    julia> model3_3_dmin,model3_3_rfr,model3_3_rmin = blackbox_model_3_3(
        exp_dmin["Data"]["dltime"],exp_dmin["Fitting"]["Attenuation"][sel_dmin,:],sel_dmin,
        exp_rfr["Data"]["dltime"],exp_rfr["Fitting"]["Attenuation"][sel_rfr,:],sel_rfr,
        exp_rmin["Data"]["dltime"],exp_rmin["Fitting"]["Attenuation"][sel_rmin,:],sel_rmin,
        coha,
        dt=dt
    )
```

