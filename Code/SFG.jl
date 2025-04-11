using Statistics, LsqFit

"""
    Fit the reference spectrum.
    Parameters:
        v             : Wavenumber vector
        spectrum      : Measured spectrum to be fitted
        initial_guess : A vector containing initial guesses for [A, omega, gamma, phi]
        lower_limit   : Lower bounds for the parameters
        upper_limit   : Upper bounds for the parameters

   Keyword Arguments:
        maxIter : Maximum number of iterations for the optimizer (default 1000)
        g_tol   : Tolerance for the gradient (default 1e-8)
        x_tol   : Tolerance for the parameters (default 1e-8)
        error   : If true, returns both the fit parameters and the full result

        julia> fit_reference_general(
            v,
            spectrum,
            initial_guess,
            lower_limit,
            upper_limit;
            maxIter = 1000,	
            g_tol = 1e-8,
            x_tol = 1e-8,
            error = false
        )
"""
function fit_reference_general(v,spectrum,initial_guess,lower_limit,upper_limit;
    maxIter = 1000,	
    g_tol = 1e-8,
    x_tol = 1e-8,
    error = false
    )

    # Determine the number of modes (n) by dividing the length of initial_guess by 4.
        n = length(initial_guess)÷4
    # Use the constant phase relation proposed by Linke. 
    # Therefore we need to identify the different phases in the initial guesses
    # When two or more phases are the same, they will be the same after fitting.
        phi = initial_guess[3n+1:4n]
        different_phi = unique(phi[phi.>0])

    # For each unique phase value, determine all indices in phi that equal that value.
        hold_phase = [findall(==(val), phi) for val in different_phi]

    # Now extract the rest of the parameters from the initial guess.
    # Amplitudes
        A =initial_guess[1:n]
    # Frequencies
        omega = initial_guess[n+1:2n]
    # Damping factors
        gamma = initial_guess[2n+1:3n]
    # Now we construct a new initial guess, where only the unique phases are kept. 
    initial_guess_01= [A;omega;gamma;different_phi]

    # The sfg_same_phase function will be used to fit the data. It has to be defined here,
    # in order to use the hold_phase variable and there keep the initial grouped phases
    function sfg_same_phase(v,p)
        # Extract the parameters from p
        A =p[1:n]
        omega = p[n+1:2n]
        gamma = p[2n+1:3n]
        _phi = p[3n+1:end]
        # Re assign the phases to the correct indices in phi_01
        # This is done by using the hold_phase variable, which contains the indices of the
        # different phases in the original phi vector.
        phi_01 = zeros(length(gamma))
        for (i, phase_indices) in enumerate(hold_phase)
            phi_01[phase_indices] .= _phi[i]
        end
        
        # Now we can calculate the signal using the parameters
        # The signal is a complex number, so we need to use complex(0) to initialize it.
        signal = fill(complex(0),length(v))
        
        # The signal is calculated by summing over all modes.
        for i in 1:n
            signal += A[i] ./(v .-omega[i] .+ 1im*gamma[i]) .*exp(1im*phi_01[i])
        end

        # The Intensity is the absolute value of the signal squared.
        return abs2.(signal)
    end

    # Now we can fit the data using the curve_fit function from LsqFit.
    # The initial guess is the initial_guess_01 vector, which contains the amplitudes,
    result = curve_fit(
        sfg_same_phase,       # model function
        v,                    # wavenumber vector
        spectrum,             # measured spectrum
        initial_guess_01,     # modified initial guess
        lower =lower_limit,   # lower bounds for the parameters
        upper =  upper_limit, # upper bounds for the parameters
        maxIter = maxIter,    # maximum number of iterations
        g_tol = g_tol,        # tolerance for the gradient
        x_tol = x_tol         #
    )

    # Since we used the constant phase relation, we need to reassign the phases to the correct indices in phi_01 again.
    phi_01 = zeros(length(phi))
    for (i, phase_indices) in enumerate(hold_phase)
        phi_01[phase_indices] .= result.param[3n + i]
    end
    # Now we can recombine the optimized parameters and the reconstructed phase vector
    fit_params = [result.param[1:3n];phi_01]

    # Post-processing: Adjust parameters that hit the bounds
    #
    # For each mode, if the parameter difference between the fitted value and its 
    # bound is zero, modify the parameter slightly (by ±1e-3) to avoid potential issues with LsqFit.

    for i in 1:n

        _A = fit_params[i]
        _A_lower = lower_limit[i]
        _A_upper = upper_limit[i]

        _omega = fit_params[i+n]
        _omega_lower = lower_limit[i+n]
        _omega_upper = upper_limit[i+n]

        _gamma = fit_params[i+2n]
        _gamma_lower = lower_limit[i+2n]
        _gamma_upper = upper_limit[i+2n]

        # Prepare arrays of parameters and their bounds for convenience
        _p = [_A,_omega,_gamma]
        _p_lower = [_A_lower,_omega_lower,_gamma_lower]
        _p_upper = [_A_upper,_omega_upper,_gamma_upper]

        # Compute the differences to the lower and upper bounds
        _d_plower = _p .- _p_lower
        _d_pupper = _p .- _p_upper

        # Adjust parameters that are exactly at the lower bound.
        if any(==(0),_d_plower)
            args = findall(==(0),_d_plower)
            for arg in args
                if arg == 1
                    fit_params[i] += 1e-3
                elseif arg == 2
                    fit_params[i+n] += 1e-3
                elseif arg == 3
                    fit_params[i+2n] += 1e-3
                end
            end
        end
        # Adjust parameters that are exactly at the upper bound.
        if any(==(0),_d_pupper)
            args = findall(==(0),_d_pupper)
            for arg in args
                if arg == 1
                    fit_params[i] -=  1e-3
                elseif arg == 2
                    fit_params[i+n] -= 1e-3
                elseif arg == 3
                    fit_params[i+2n] -=  1e-3
                end
            end
        end  
    end     

    # Optionally return both the fit parameters and the complete result structure.
    if error 
        return fit_params, result
    end

    return fit_params
end

"""
   Fit the signal and difference spectra with constant reference spectrum. This function uses the 
   reference fit parameters (from fit_reference_general) and further refines the 
   signal fit by taking into account an attenuation factor that is adjusted over a series 
   of measurements. Therefore we assume that in the signal spectra all parameters are kept konstant except
   the amplitude A. The function uses the LsqFit package to perform the fitting.

   Parameters:
        initial_guess : Initial guesses for the model parameters
        data          : The data dictionary

   Keyword Arguments:
        hold        : Vector indicating which modes are held constant (1) or fitted (0)
        lower_A     : Lower bound for amplitude A
        lower_omega : Lower bound for omega (offset from initial guess by -5.0)
        lower_gamma : Lower bound for gamma
        lower_phi   : Lower bound for phase (calculated from unique phases > 0)
        upper_A     : Upper bound for amplitude A
        upper_omega : Upper bound for omega (offset from initial guess by +5.0)
        upper_gamma : Upper bound for gamma
        upper_phi   : Upper bound for phase (set to 2π)

    Returns:
        A tuple containing the fit parameters from the reference fit and an array of fit
        results from subsequent attenuated signal fits.
    julia> fit_signal_const_ref_attenuation(initial_guess,data;
        hold=zeros(length(initial_guess)÷4),
        lower_A = fill(0.1,length(initial_guess)÷4),
        lower_omega = initial_guess[((length(initial_guess)÷4)+1):(2*(length(initial_guess)÷4))] .-5.0,
        lower_gamma = fill(5.5,length(initial_guess)÷4),
        lower_phi  = zeros(length(unique(initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))][map(x-> x>0,initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))])]))),
        upper_A     = fill(20,length(initial_guess)÷4),
        upper_omega = initial_guess[((length(initial_guess)÷4)+1):(2*(length(initial_guess)÷4))] .+5.0,
        upper_gamma = fill(50,length(initial_guess)÷4),
        upper_phi  = fill(2*pi,length(unique(initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))][map(x-> x>0,initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))])]))),
        )

"""
function fit_signal_const_ref_attenuation(initial_guess,data;
    hold=zeros(length(initial_guess)÷4),         
    lower_A = fill(0.1,length(initial_guess)÷4),
    lower_omega = initial_guess[((length(initial_guess)÷4)+1):(2*(length(initial_guess)÷4))] .-5.0,
    lower_gamma = fill(5.5,length(initial_guess)÷4),
    lower_phi  = zeros(length(unique(initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))][map(x-> x>0,initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))])]))),
    upper_A     = fill(20,length(initial_guess)÷4),
    upper_omega = initial_guess[((length(initial_guess)÷4)+1):(2*(length(initial_guess)÷4))] .+5.0,
    upper_gamma = fill(50,length(initial_guess)÷4),
    upper_phi  = fill(2*pi,length(unique(initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))][map(x-> x>0,initial_guess[(3(length(initial_guess)÷4)+1):(4*(length(initial_guess)÷4))])]))),
    )

    # Construct the overall lower and upper bounds by concatenating the parameter bounds.
        lower_limit = [lower_A;lower_omega;lower_gamma;lower_phi] 
        upper_limit = [upper_A;upper_omega;upper_gamma;upper_phi]

    # Extract the wavenumber vector
        v = data["wavenumber"]
    # We want to create ROI where we want the fit to be more precise. We wil
    # First we initialise the weight (wt) vector with ones. The length
    # of wt is twice the length of v since we fit the signal spectrum and the difference spectrum
    # (signal - reference) per timestep.
        wt = ones(2*length(v))
    # In the spectral region between 2810 and 3000 cm-1 we want to give more weight to the fit.
        roi = (2810,3000)
        roi_arg= findall(x-> roi[1]<=x<=roi[2],v)
    # Weight of 2 in the ROI for the signal spectrum
        wt[roi_arg] .= 2.0
    # Weight of 3 in the ROI for the difference spectrum (signal - reference)
        wt[roi_arg .+ length(v)] .= 3.0

    # Now we can fit the reference spectrum using the fit_reference_general function.
        fit_params = fit_reference_general(v, mean(data["ref matrix"],dims=1)[:] ,initial_guess,lower_limit,upper_limit)
    # We concatenate the signalspectrum at time i with the difference spectrum (signal - reference) at time i.
        m_sigdiff = [[data["sig matrix"][i,:];data["sig matrix"][i,:] - mean(data["ref matrix"],dims=1)[:]] for i in 1:size(data["sig matrix"],1)]
    # Count the number of modes to adjust (when hold value is 0).
        fit_n = count(==(0),hold)
    # Initialize an empty array to collect fit results for each signal.
        fits = []
    # Model function to calculate the signal and difference spectrum by only varying the Amplitude quantified by
    # attenuation coefficent a. The rest of the parameters are kept constant.
    function sfg_attenuation(v,p::Vector{T}) where T<:Number

        p0 = fit_params                     # Reference fit parameters
        a = ones(T,n)                       # Atennuation coefficients set to 1.0 on default
            a[hold .== 0] .= p[1:fit_n]     # Set a[i] to p[i] if hold[i] == 0
        A = p0[1:n]                         # Amplitudes from reference fit
        omega = p0[n+1:2n]                  # Frequencies from reference fit
        gamma = p0[2n+1:3n]                 # Damping factors from reference fit
        phi  = p0[3n+1:4n]                  # Phases from reference fit
        ref  = fill(complex(0),length(v))   # Initialize reference spectrum
        signal = fill(complex(0),length(v)) # Initialize signal spectrum

        # Calculate the reference spectrum and the signal spectrum
        # The spectra are calculated by summing over n modes.
            for i in 1:n
                ref    +=         A[i]  ./(v .-omega[i] .+ 1im*gamma[i]) .*exp(1im*phi[i])
                signal += (a[i] * A[i]) ./(v .-omega[i] .+ 1im*gamma[i]) .*exp(1im*phi[i])
            end
        # The Intensity is the absolute value of the signal squared.
            ref = abs2.(ref)
            signal = abs2.(signal)
        # The difference spectrum is calculated by subtracting the reference spectrum from the signal spectrum.
            diff = signal - ref
        return [signal;diff]
    end

    # Now we can fit the signal and difference spectra using the sfg_attenuation function
    # We loop over all timesteps
    for i in 1:size(data["ref matrix"],1)
        # For the first iteration, we set the lower and upper bounds for the attenuation coefficients to 0.95 and 1.05.
        # Since we expect no big change  for the first iteration
        # The [hold .== 0] part is used to select the parameters that are not held constant.
        if  i == 1 
            lower_a = fill(0.95,length(fit_params[1:n][hold .== 0])) # lower bound for attenuation
            upper_a = fill(1.05,length(fit_params[1:n][hold .== 0])) # upper bound for attenuation
            lower= lower_a
            upper= upper_a
            # Use the curve_fit function from LsqFit Package
            fit = curve_fit(
                sfg_attenuation,                               # model function
                data["wavenumber"],                            # wavenumber vector                                    
                m_sigdiff[i],                                  # signal and difference spectrum at time i                         
                wt,                                            # weight vector
                fill(1.0,length(fit_params[1:n][hold .== 0])), # initial guess for attenuation 
                lower=lower,                                   # lower bounds for the parameters
                upper=upper                                    # upper bounds for the parameters
            )
            # Push the fit result into the fits array
            push!(fits,fit)
        else
            # For the subsequent iterations, we set the lower and upper bounds for the attenuation dynamically
            # based on the previous fit parameters. The lower bounds are set to the fitted values minus 0.5 (d_a), 
            # and the upper bounds are set to be no bigger than 3.0 (threshold_a)
                d_a = 0.5
                threshold_a = 3.0
            # Calculate the dynamical lower and upper bounds for the attenuation. Check if the lower bound is less 
            # than 0.0 and set it to 0.0 if so or if the upper bound is greater than threshold_a then set it to threshold_a.
                lower_a =[a-d_a < 0 ? 0.0 : a-d_a for a in fits[i-1].param[1:fit_n]]
                upper_a = [a+d_a > threshold_a ? threshold_a : a+d_a for a in fits[i-1].param[1:fit_n]]
                lower= lower_a
                upper= upper_a
            # Use the curve_fit function from LsqFit Package
            fit = curve_fit(
                sfg_attenuation,                                # model function
                data["wavenumber"],                             # wavenumber vector                       
                m_sigdiff[i],                                   # signal and difference spectrum at time i 
                wt,                                             # weight vector
                fits[i-1].param,                                # initial guess based on previous fit
                lower=lower,                                    # lower bounds for the parameters
                upper=upper                                     # upper bounds for the parameters
            )
            # Push the fit result into the fits array
            push!(fits,fit)
        end
    end
    return fit_params,fits
end