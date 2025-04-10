

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
    time=10.0
    )
    
    lower_limit = [lower_A;lower_omega;lower_gamma;lower_phi] 

    upper_limit = [upper_A;upper_omega;upper_gamma;upper_phi]

    v = data["wavenumber"]
    wt = ones(2*length(v))
    roi = (2810,3000)
    roi_arg= findall(x-> roi[1]<=x<=roi[2],v)
    wt[roi_arg] .= 2.0
    wt[roi_arg .+ length(v)] .= 3.0

    fit_params = fit_reference_general(v, mean(data["ref matrix"],dims=1)[:] ,initial_guess,lower_limit,upper_limit)
    m_sigdiff = [[data["sig matrix"][i,:];data["sig matrix"][i,:] - mean(data["ref matrix"],dims=1)[:]] for i in 1:size(data["sig matrix"],1)]

    fit_n = count(==(0),hold)
    fits = []
    function sfg_attenuation(v,p::Vector{T}) where T<:Number
        p0 = fit_params 
        a = ones(T,n)
            a[hold .== 0] .= p[1:fit_n]
        A = p0[1:n]
        omega = p0[n+1:2n]
        gamma = p0[2n+1:3n]
        phi  = p0[3n+1:4n]
        ref  = fill(complex(0),length(v))
        signal = fill(complex(0),length(v))

        for i in 1:n
            ref    +=         A[i]  ./(v .-omega[i] .+ 1im*gamma[i]) .*exp(1im*phi[i])
            signal += (a[i] * A[i]) ./(v .-omega[i] .+ 1im*gamma[i]) .*exp(1im*phi[i])
        end
        ref = abs2.(ref)
        signal = abs2.(signal)
        diff = signal - ref
        return [signal;diff]
    end

    for i in 1:size(data["ref matrix"],1)

        if  i == 1 
            lower_a = fill(0.95,length(fit_params[1:n][hold .== 0]))
            upper_a = fill(1.05,length(fit_params[1:n][hold .== 0]))
            lower= lower_a
            upper= upper_a
            fit = curve_fit(
                sfg_attenuation,
                data["wavenumber"],
                m_sigdiff[i],
                wt,
                fill(1.0,length(fit_params[1:n][hold .== 0])),
                lower=lower,
                upper=upper
            )
            push!(fits,fit)
        else
            d_a = 0.5
            threshold_a = 3.0
            lower_a =[a-d_a < 0 ? 0.0 : a-d_a for a in fits[i-1].param[1:fit_n]]
            upper_a = [a+d_a > threshold_a ? threshold_a : a+d_a for a in fits[i-1].param[1:fit_n]]
            lower= lower_a
            upper= upper_a
            fit = curve_fit(
                sfg_attenuation,
                data["wavenumber"],
                m_sigdiff[i],
                wt,
                fits[i-1].param,
                lower=lower,
                upper=upper
            )
            push!(fits,fit)
        end
    end
    return fit_params,fits
end