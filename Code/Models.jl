using DataFrames
include("Reservoir Model Struct.jl")
"""
    Model 3.3
    Model 3.3 is an extension of Model 3.2. Certain rates appear to be sufficiently slow that they can be effectively removed from the simulation without 
    significantly affecting accuracy.
    The parameters that have been removed or fixed in this model are listed explicitly under model.Dictionary.
    Modell 3.2
    From wavelength-scan experiments, it became clear that pumping at the wavelength labeled "rmin" actually excites two modes simultaneously 
    (rmina and rminb). Therefore, an additional pump pulse parameter has been introduced. This parameter is set to zero whenever "rmin" is not 
    actively pumped.
    Model 3.2 is considered an extension and refinement of Model 3.1.
    Modell 3.1:
    Model 3.1 is derived from Model 2.3. The primary shortcoming of Model 2.3 was its inability to accurately reproduce the fast dynamics 
    observed in the "rfr" mode. Model 3.1 addresses this issue by allowing forward and backward rates involving the "rfr" mode to differ (asymmetric transitions).
    As in Model 2.3, parameters p[1:2] describe the laser excitation processes. Parameters p[3:9] represent the relaxation rates of the modes, 
    whereas p[10:50] specify the coupling rates between different modes. Finally, parameters p[51:57] define the degeneracy of each mode.
        julia> model3_3(
                t,                  # model time vector
                p::Vector{T};       # model parameters
                mode = :dmin,       # which mode is pumped? :dmin, :rmin, or :rfr
                dt = 0.1            # time step size for the simulation
               ) where T<:Number

"""
function model3_3(
    t,
    p::Vector{T};
    mode = :dmin,
    dt = 0.1
) where T<:Number

    #---------------------------------------------------------
    # 1) We have 7 excited states. We'll store their populations
    #    in Popu[:,1] through Popu[:,7].
    #    We'll also store a "bleach" signal in Bleach[:,1..7].
    #---------------------------------------------------------
    n = 7  # number of states
    Popu   = zeros(T, length(t), n)
    Bleach = zeros(T, length(t), n)
    control = zeros(T, length(t), 4)

    #---------------------------------------------------------
    # 2) Unpack parameters from p:
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    α_rb     = 0      # amplitude scaling of laser pulse for rminb/u
    if mode == :rmin
        α_rb = p[2]
    end
    t_0      = p[3]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width to observed pulse width of pump pulse

    # Relaxation rates to ground state 
    relax_rplus       = 1/p[4]
    #relax_dminusomega = 1/p[5]
    #relax_dfr         = 1/p[6]
    relax_dminus      = 1/p[5]
    #relax_rfr         = 1/p[4]
    relax_rminusb     = 1/p[6]
    relax_rminusa     = 1/p[7]

    # Coupling constants between states 
    # For simplicity, I'll rename them 'k_x_y' to emphasize
    # they are coupling rates. Each will get multiplied by
    # the product of degeneracies g_x*g_y for symmetrical flow.

    k_rplus_dminusomega   = 1/p[8]
    k_rplus_dfr           = 1/p[9]
    #k_rplus_dminus        = 1/p[13]
    #k_rplus_rfr           = 1/p[10]
    k_rplus_rminusb       = 1/p[10]
    #k_rplus_rminusa       = 1/p[16]

    k_dminusomega_rplus   = k_rplus_dminusomega
    #k_dminusomega_dfr     = 1/p[10]
    #k_dminusomega_dminus  = 1/p[18]
    k_dminusomega_rfr     = 1/p[11]
    #k_dminusomega_rminusb = 1/p[20]
    k_dminusomega_rminusa = 1/p[12]

    k_dfr_rplus           = k_rplus_dfr
    #k_dfr_dminusomega     = k_dminusomega_dfr
    #k_dfr_dminus          = 1/p[22]
    #k_dfr_rfr             = 1/p[23]
    k_dfr_rminusb         = 1/p[13]
    #k_dfr_rminusa         = 1/p[15]

    #k_dminus_rplus        = k_rplus_dminus
    #k_dminus_dminusomega  = k_dminusomega_dminus
    #k_dminus_dfr          = k_dfr_dminus
    k_dminus_rfr          = 1/p[14]
    #k_dminus_rminusb      = 1/p[27]
    #k_dminus_rminusa      = 1/p[28]

    #k_rfr_rplus           = 1/p[17]
    k_rfr_dminusomega     = 1/p[15]
    #k_rfr_dfr             = 1/p[31]
    k_rfr_dminus          = 1/p[16]
    k_rfr_rminusb         = 1/p[17]
    k_rfr_rminusa         = 1/p[18]

    k_rminusb_rplus       = k_rplus_rminusb
    #k_rminusb_dminusomega = k_dminusomega_rminusb
    k_rminusb_dfr         = k_dfr_rminusb
    #k_rminusb_dminus      = k_dminus_rminusb
    k_rminusb_rfr         = 1/p[19]
    #k_rminusb_rminusa     = 1/p[19]

    #k_rminusa_rplus       = k_rplus_rminusa
    k_rminusa_dminusomega = k_dminusomega_rminusa
    #k_rminusa_dfr         = k_dfr_rminusa
    #k_rminusa_dminus      = k_dminus_rminusa
    k_rminusa_rfr         = 1/p[20]
    #k_rminusa_rminusb     = k_rminusb_rminusa

    # 3) Degeneracy of the states:
    # Fixed values after Model 3.2
        g_rplus   = 1 #p[38]
        g_d_omega = 2
        g_dfr     = 12
        g_dminus  = 16
        g_rfr     = 1
        g_rminb   = 3
        g_rmina   = 2 #p[44]

    # Store them in an array:
    #   1->rplus, 2->dminusomega, 3->dfr, 4->dminus, 
    #   5->rfr,   6->rminusb,     7->rminusa
        G = [g_rplus, g_d_omega, g_dfr, g_dminus, g_rfr, g_rminb, g_rmina]

    #---------------------------------------------------------
    # 4) Precompute the pulse shape at each time point
        f_g    = gauss.(t,Ref([σ,t_0]))

    #---------------------------------------------------------
    # 5) Mode-specific initial conditions:
    #    :dmin, :rmin, or :rfr  -> which state is initially pumped?
    #---------------------------------------------------------
        pump_first = α * f_g[1] * dt

    # Initially, one state is pumped:
        Popu[1,1] = (mode == :rplus)     ? 1/G[1] * pump_first         : 0
        Popu[1,2] = (mode == :dminomega) ? 1/G[2] * pump_first         : 0
        Popu[1,3] = (mode == :dfr)       ? 1/G[3] * pump_first         : 0
        Popu[1,4] = (mode == :dmin)      ? 1/G[4] * pump_first         : 0
        Popu[1,5] = (mode == :rfr)       ? 1/G[5] * pump_first         : 0
        Popu[1,6] = (mode == :rmin)      ? 1/G[6] * α_rb * f_g[1] * dt : 0
        Popu[1,7] = (mode == :rmin)      ? 1/G[7] * pump_first         : 0
        control[1,1] = pump_first + α_rb * f_g[1] * dt 
        control[1,2] = G[1] * Popu[1,1] + (
            + G[2] * Popu[1,2] 
            + G[3] * Popu[1,3] 
            + G[4] * Popu[1,4] 
            + G[5] * Popu[1,5] 
            + G[6] * Popu[1,6] 
            + G[7] * Popu[1,7]
        )


    #---------------------------------------------------------
    # 6) Time stepping with explicit Euler method:
    #    We'll update all 7 states at each time step.
    #
    #    Key change for symmetrical coupling:
    #      * For a coupling k_xy between state x and y,
    #        we multiply by (G[x]*G[y]) in BOTH equations:
    #
    #        Popu[i,x] += + k_xy * G[x]*G[y]*(Popu[i-1,y] - Popu[i-1,x])
    #        Popu[i,y] += + k_xy * G[x]*G[y]*(Popu[i-1,x] - Popu[i-1,y])
    #
    #    The same factor appears in both directions.
    #---------------------------------------------------------
    for i in 2:length(t)

        # Precompute the laser amplitude at this timestep
        laser_here = α * f_g[i]

        #-------------
        # State 1: rplus 
        #-------------
        # If this is the pumped state in :rplus mode, we add laser_here
        # (same code runs for all modes, but laser_here=0 if not used)
        pump_rplus = (mode == :rplus) ?  laser_here : 0
        Popu[i,1] = Popu[i-1,1] + 1/G[1] * (

            # Pumping if this mode is :rplus
            + pump_rplus

            # Relaxation to ground
            - G[1] * relax_rplus * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rplus_dminusomega * (G[1]*G[2]) * Popu[i-1,1]
            + k_dminusomega_rplus * (G[2]*G[1]) * Popu[i-1,2] 

            # Coupling with dfr (state 3)
            - k_rplus_dfr         * (G[1]*G[3]) * Popu[i-1,1]
            + k_dfr_rplus         * (G[3]*G[1]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            #- k_rplus_dminus      * (G[1]*G[4]) * Popu[i-1,1]
            #+ k_dminus_rplus      * (G[4]*G[1]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            #- k_rplus_rfr         * (G[1]*G[5]) * Popu[i-1,1]
            #+ k_rfr_rplus         * (G[5]*G[1]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rplus_rminusb     * (G[1]*G[6]) * Popu[i-1,1]
            + k_rminusb_rplus     * (G[6]*G[1]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            #- k_rplus_rminusa     * (G[1]*G[7]) * Popu[i-1,1]
            #+ k_rminusa_rplus     * (G[7]*G[1]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 2: dminusomega
        #-------------
        pump_dminomega = (mode == :dminomega) ? laser_here : 0
        Popu[i,2] = Popu[i-1,2] + 1/G[2] * (

            # Pumping if this mode is :dminomega
            + pump_dminomega

            # Relaxation
            #- G[2] * relax_dminusomega * Popu[i-1,2]

            # Coupling with rplus (state 1)
            - k_dminusomega_rplus   * (G[2]*G[1]) * Popu[i-1,2]
            + k_rplus_dminusomega   * (G[1]*G[2]) * Popu[i-1,1]

            # Coupling with dfr (state 3)
            #- k_dminusomega_dfr     * (G[2]*G[3]) * Popu[i-1,2]
            #+ k_dfr_dminusomega     * (G[3]*G[2]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            #- k_dminusomega_dminus  * (G[2]*G[4]) * Popu[i-1,2]
            #+ k_dminus_dminusomega  * (G[4]*G[2]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dminusomega_rfr     * (G[2]*G[5]) * Popu[i-1,2]
            + k_rfr_dminusomega     * (G[5]*G[2]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            #- k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]
            #+ k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]
            + k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 3: dfr
        #-------------
        pump_dfr = (mode == :dfr) ? laser_here : 0
        Popu[i,3] = Popu[i-1,3] + 1/G[3] * (

            # Pumping if this mode is :dfr
            + pump_dfr

            # Relaxation
            #- G[3] * relax_dfr * Popu[i-1,3]

            # Coupling with rplus (state 1)
            - k_dfr_rplus       * (G[3]*G[1]) * Popu[i-1,3]
            + k_rplus_dfr       * (G[1]*G[3]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            #- k_dfr_dminusomega * (G[3]*G[2]) * Popu[i-1,3]
            #+ k_dminusomega_dfr * (G[2]*G[3]) * Popu[i-1,2]

            # Coupling with dminus (state 4)
            #- k_dfr_dminus      * (G[3]*G[4]) * Popu[i-1,3]
            #+ k_dminus_dfr      * (G[4]*G[3]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            #- k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]
            #+ k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dfr_rminusb     * (G[3]*G[6]) * Popu[i-1,3]
            + k_rminusb_dfr     * (G[6]*G[3]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            #- k_dfr_rminusa     * (G[3]*G[7]) * Popu[i-1,3]
            #+ k_rminusa_dfr     * (G[7]*G[3]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 4: dminus
        #-------------
        pump_dmin = (mode == :dmin) ? laser_here : 0
        Popu[i,4] = Popu[i-1,4] + 1/G[4] * (

            # Pumping if this mode is :dmin
            + pump_dmin

            # Relaxation
            - G[4] * relax_dminus * Popu[i-1,4]

            # Coupling with rplus (state 1)
            #- k_dminus_rplus       * (G[4]*G[1]) * Popu[i-1,4]
            #+ k_rplus_dminus       * (G[1]*G[4]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            #- k_dminus_dminusomega * (G[4]*G[2]) * Popu[i-1,4]
            #+ k_dminusomega_dminus * (G[2]*G[4]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            #- k_dminus_dfr         * (G[4]*G[3]) * Popu[i-1,4]
            #+ k_dfr_dminus         * (G[3]*G[4]) * Popu[i-1,3]

            # Coupling with rfr (state 5)
            - k_dminus_rfr         * (G[4]*G[5]) * Popu[i-1,4]
            + k_rfr_dminus         * (G[5]*G[4]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            #- k_dminus_rminusb     * (G[4]*G[6]) * Popu[i-1,4]
            #+ k_rminusb_dminus     * (G[6]*G[4]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            #- k_dminus_rminusa     * (G[4]*G[7]) * Popu[i-1,4]
            #+ k_rminusa_dminus     * (G[7]*G[4]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 5: rfr
        #-------------
        pump_rfr = (mode == :rfr) ? laser_here : 0
        Popu[i,5] = Popu[i-1,5] + 1/G[5] * (
            # Pumping if this mode is :rfr
            + pump_rfr
            
            # Relaxation
            #- G[5] * relax_rfr * Popu[i-1,5]

            # Coupling with rplus (state 1)
            #- k_rfr_rplus       * (G[5]*G[1]) * Popu[i-1,5]
            #+ k_rplus_rfr       * (G[1]*G[5]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rfr_dminusomega * (G[5]*G[2]) * Popu[i-1,5]
            + k_dminusomega_rfr * (G[2]*G[5]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            #- k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]
            #+ k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rfr_dminus      * (G[5]*G[4]) * Popu[i-1,5]
            + k_dminus_rfr      * (G[4]*G[5]) * Popu[i-1,4]

            # Coupling with rminusb (state 6)
            - k_rfr_rminusb     * (G[5]*G[6]) * Popu[i-1,5]
            + k_rminusb_rfr     * (G[6]*G[5]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rfr_rminusa     * (G[5]*G[7]) * Popu[i-1,5]
            + k_rminusa_rfr     * (G[7]*G[5]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 6: rminusb or u 
        #-------------
        #pump_rminb = (mode == :rminb) ? laser_here : 0
        pump_rminb  = α_rb * f_g[i] * dt
        Popu[i,6] = Popu[i-1,6] + 1/G[6] * (

            # Pumping if this mode is :rminb
            + pump_rminb

            # Relaxation
            - G[6] * relax_rminusb * Popu[i-1,6]

            # Coupling with rplus (state 1)
            - k_rminusb_rplus       * (G[6]*G[1]) * Popu[i-1,6]
            + k_rplus_rminusb       * (G[1]*G[6]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            #- k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]
            #+ k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusb_dfr         * (G[6]*G[3]) * Popu[i-1,6]
            + k_dfr_rminusb         * (G[3]*G[6]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            #- k_rminusb_dminus      * (G[6]*G[4]) * Popu[i-1,6]
            #+ k_dminus_rminusb      * (G[4]*G[6]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusb_rfr         * (G[6]*G[5]) * Popu[i-1,6]
            + k_rfr_rminusb         * (G[5]*G[6]) * Popu[i-1,5]

            # Coupling with rminusa (state 7)
            #- k_rminusb_rminusa     * (G[6]*G[7]) * Popu[i-1,6]
            #+ k_rminusa_rminusb     * (G[7]*G[6]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 7: rminusa
        #-------------
        pump_rmin = (mode == :rmin) ? laser_here : 0
        Popu[i,7] = Popu[i-1,7] + 1/G[7] * (

            # Pumping if this mode is :rmin
            + pump_rmin

            # Relaxation
            - G[7] * relax_rminusa * Popu[i-1,7]

            # Coupling with rplus (state 1)
            #- k_rminusa_rplus       * (G[7]*G[1]) * Popu[i-1,7]
            #+ k_rplus_rminusa       * (G[1]*G[7]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
            + k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            #- k_rminusa_dfr         * (G[7]*G[3]) * Popu[i-1,7]
            #+ k_dfr_rminusa         * (G[3]*G[7]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            #- k_rminusa_dminus      * (G[7]*G[4]) * Popu[i-1,7]
            #+ k_dminus_rminusa      * (G[4]*G[7]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusa_rfr        * (G[7]*G[5]) * Popu[i-1,7]
            + k_rfr_rminusa        * (G[5]*G[7]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            #- k_rminusa_rminusb    * (G[7]*G[6]) * Popu[i-1,7]
            #+ k_rminusb_rminusa    * (G[6]*G[7]) * Popu[i-1,6]
        ) * dt

        #-------------
        # 7) Control
        # Checks if the balance of the populations is correct at each timestep
        #-------------
        
        # 1. Integral of Laser Input 
        control[i,1] = control[i-1,1] + (laser_here + α_rb * f_g[i]) * dt
        # 2. Overall Population of all substates
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )
        # 3. Relaxed Population of all substates
        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            #+ G[2] * relax_dminusomega * Popu[i-1,2]
            #+ G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            #+ G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        # 4. Sum of Laser Input, Overall Population and Relaxed Population should equal to 0 at every timestep
        control[i,4] = control[i,1] - control[i,2] - control[i,3]
    end

    #---------------------------------------------------------
    # 8) Finally, compute the "Bleach" signal 
    #    for each state: (1 - 2*g[i]*Popu[i])^2
    #---------------------------------------------------------
    for k in 1:n
        Bleach[:,k] = (1 .- 2 .* G[k] .* Popu[:,k]) .^2
    end

    # ------------------------------------------------------
    # 9) Dictionary (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1],alpha_rb = p[2], t_0 = p[3])

    # B) Relaxation
    modes_relax = [
        "rplus", 
        #"dminusomega", 
        #"dfr", 
        "dminus", 
        #"rfr", 
        "rminusb", 
        "rminusa"
        ]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[4:7]
    )

    # C) Degeneracy
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G#p[38:44]
    )

    # D) Couplings between states
    #    We differ between forward and backward couplings
    #    Syntax: (x, y, fwd_idx, bwd_idx)
    #      => t_x_y   => p[fwd_idx]
    #      => t_y_x   => p[bwd_idx]
    couplings = [
        (:rplus,       :dminusomega, 8, 8),
        (:rplus,       :dfr,         9, 9),
        (:rplus,       :rminusb,     10, 10), # asym

        (:dminusomega, :rfr,         11, 15),  # asym
        (:dminusomega, :rminusa,     12, 12),

        (:dfr,         :rminusb,     13, 13),

        (:dminus,      :rfr,         14, 16),  # asym

        (:rfr,         :rminusb,     17, 19),  # asym
        (:rfr,         :rminusa,     18, 20),  # asym
    ]

    rows_all = Vector{Dict{Symbol,Any}}()
    for (x, y, fwd_idx, bwd_idx) in couplings
        push!(rows_all, Dict(
            :param_name => "t_$(x)_$(y)",
            :value      => p[fwd_idx],
            :ratio     => p[fwd_idx]/p[bwd_idx]

        ))
        push!(rows_all, Dict(
            :param_name => "t_$(y)_$(x)",
            :value      => p[bwd_idx],
            :ratio     => p[bwd_idx]/p[fwd_idx]
        ))
    end
    df_koppl_all = DataFrame(rows_all)
    # Function to create DataFrame for each mode
    function df_koppl_for_mode(m::Symbol)
        subrows = Vector{Dict{Symbol,Any}}()
        for (x,y, fwd_idx, bwd_idx) in couplings
            if x == m || y == m
                # x->y
                push!(subrows, Dict(
                    :param_name => "t_$(x)_$(y)",
                    :value      => p[fwd_idx],
                    :ratio     => p[fwd_idx]/p[bwd_idx]
                ))
                # y->x
                push!(subrows, Dict(
                    :param_name => "t_$(y)_$(x)",
                    :value      => p[bwd_idx],
                    :ratio     => p[bwd_idx]/p[fwd_idx]
                ))
            end
        end
        return DataFrame(subrows)
    end
    # Create DataFrame for each mode
    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end
    # E) Removed parameters
      df_entfernt = DataFrame(
        Parameter = [
            "relax_dminusomega",
            "relax_dfr",
            "relax_rfr",
            "k_rplus_dminus",
            "k_rplus_rfr",
            "k_rplus_rminusa",
            "k_dminusomega_dfr",
            "k_dminusomega_dminus",
            "k_dminusomega_rminusb",
            "k_dfr_dminus",
            "k_dfr_rfr",
            "k_dfr_rminusa",
            "k_dminus_rminusb",
            "k_dminus_rminusa",
            "k_rfr_rplus",
            "k_rfr_dfr",
            "k_rminusb_rminusa",
            "g_rplus",
            "g_d_omega",
            "g_dfr",
            "g_dminus",
            "g_rfr",
            "g_rminusb",
            "g_rminusa",
            ],
        Wert = vcat(
            fill("Entfernt",3),
            fill("Entfernt",14),
            1,
            2,
            12,
            16,
            1,
            3,
            2
            )
      )
    # Full output dictionary
    outputDict = Dict(
        :Model       => "Model 3.3",
        :fitted_mode => mode,
        :Parameter   => Dict(
            :All        => p,
            :Puls       => df_puls,
            :Relaxation => df_relax,
            :Entartung  => df_entartung,
            :Kopplung   => dict_koppl,
            :Entfernt   => df_entfernt
        ),
        :Zeitschritt => dt
    )

    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end

"""
    The following function can fit the model 3.3 to the data. For more information on the model, see the model function model3_3.
    The function uses a global optimizer of the package "BlackBoxOptim.jl". It is recommend to use the :adaptive_de_rand_1_bin_radiuslimited algorithm 
    to have a robustfree fit. 
    This function fits three data sets to one set of parameters. The data sets are:
        - time_dmin, traces_dmin, use_dmin: data set for the dmin pumped experiment
        - time_rfr, traces_rfr, use_rfr: data set for the rfr pumped experiment
        - time_rmin, traces_rmin, use_rmin: data set for the rmin pumped experiment
    The time argument is the time vector of the corresponding data set, the traces argument is the data set itself 
    of the used traces and the use argument is a vector which contains the indices of the traces which should be used for the fit.
    The coha argument is an array. coha[1] contains the index of the mode where we observed the coherent artifact in the rfr pumped experiment.
    coha[2] contains the indices of the time_vector where the coherent artifact is observed.
    The function returns three ReservoirModel objects for the three data sets.
        julia> blackbox_model_3_3(
                time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
                time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
                time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
                coha;                           # cohrent artifact array
                dt=0.1                          # time step
               )
"""
function blackbox_model_3_3(
    time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
    time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
    time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
    coha;                           # cohrent artifact array
    dt=0.1                          # time step
)
    # 1) Define the search space
    lb = vcat(
        [0.0, 0.0, -20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [0.0, 0.0, -20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [0.0, 0.0, -20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(1.0,4),  # relaxation_params   
        fill(1.0,13), # couplings_params
    )
    ub = vcat(
        [50.0, 1e-6, 20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [50.0, 1e-6, 20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [50.0, 50.0, 20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(5000.0,4),  # relaxation_params   
        fill(5000.0,13), # couplings_params
    )
    # 2) Construct time vectors for each experiment
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    # 3) Define the cost function
    # The function "cost" computes a total cost based on the weighted squared differences
    # between the model predictions and experimental traces. It includes penalties to
    # constrain the simulated populations.

    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1]; α_rb_dmin  = p[2];  t0_dmin  = p[3]
        α_rfr   = p[4]; α_rb_rfr   = p[5];  t0_rfr   = p[6]
        α_rmin  = p[7]; α_rb_rmin  = p[8];  t0_rmin  = p[9]

        relax_p = p[10:13]     # 4 relaxation
        koppl_p = p[14:26]     # 13 couplings
        g     = 1#p[28:end]    # degeneracy is for this model hardcoded in the model function

        # 2) Construct the models
        # With the use argument we select a subset of the model to fit to the data

        # dmin pumped experiment
        model_dmin = model3_3(
            t_model_dmin,
            vcat(α_dmin,α_rb_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # dmin subset
        model_dmin_used = model_dmin.Bleach[:, use_dmin]
        # rfr pumped experiment
        model_rfr = model3_3(
            t_model_rfr,
            vcat(α_rfr, α_rb_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # rfr subset
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        # rmin pumped experiment
        model_rmin = model3_3(
            t_model_rmin,
            vcat(α_rmin,α_rb_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # rmin subset
        model_rmin_used = model_rmin.Bleach[:, use_rmin]

        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
        #   Also, replace NaN and Inf values with a large number (1e12)
        pred_dmin = hcat(v_pred_dmin...)
        for i in 1:length(pred_dmin)
            if isinf(pred_dmin[i]) || isnan(pred_dmin[i])
                pred_dmin[i] = 1e12
            end
        end
        pred_rfr  = hcat(v_pred_rfr...)
        for i in 1:length(pred_rfr)
            if isinf(pred_rfr[i]) || isnan(pred_rfr[i])
                pred_rfr[i] = 1e12
            end
        end
        pred_rmin = hcat(v_pred_rmin...)
        for i in 1:length(pred_rmin)
            if isinf(pred_rmin[i]) || isnan(pred_rmin[i])
                pred_rmin[i] = 1e12
            end
        end
        # 5) Select time steps without coherence artifact in rfr pumped experiment
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights 
        # Use the time_weight function to create special weights for short delays, since we observe 
        # fast dynamics in the beginning of the experiment.
        w_dmin = ones(size(pred_dmin))
        w_pertrace_dmin = [
            5.0,  # rplus
            5.0,  # rfr
            2.0,  # rminusb
            6.0   # rminusa

        ]
        for m in 1:size(w_dmin,1)
            for n in 1:size(w_dmin,2)
                w_dmin[m,n] = time_weight(time_dmin[m],weight=3) * w_pertrace_dmin[n]
            end
        end
    
        w_rfr = ones(size(pred_rfr))
        w_pertrace_rfr = [
            5.0,  # rplus
            1.0,  # dfr
            1.0,  # rfr
            5.0   # rminusa
        ]
        for m in 1:size(w_rfr,1)
            for n in 1:size(w_rfr,2)
                if (n == coha[1]) && (m ∉ selected_idx_rfr)
                    w_rfr[m,n] = 0.0
                else
                    w_rfr[m,n] = time_weight(time_rfr[m],weight=3) * w_pertrace_rfr[n]
                end
            end
        end
        w_rmin = ones(size(pred_rmin))
        w_pertrace_rmin = [
            5.0,  # rplus
            1.0,  # dminusomega
            7.0,  # rfr
            2.0,  # rminusb
            5.0   # rminusa
        ]
        for m in 1:size(w_rmin,1)
            for n in 1:size(w_rmin,2)
                w_rmin[m,n] = time_weight(time_rmin[m],weight=3) * w_pertrace_rmin[n]
            end
        end
        
        # 7) Since we sum over all Data points, we need to normalize the weights
        # to the length of the given dataset
            length_dmin = length(time_dmin)
            length_rfr = length(time_rfr)
            length_rmin = length(time_rmin)
            norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
            norm_weight ./= maximum(norm_weight)
            cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
            cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
            cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 8) penalties
        # We define a maximum Population with max_P
        # If the Population is above this value, we add a penalty to the cost function
            max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
            P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
                norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
                norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
            total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty

        return total_cost
    end

    # Use bboptimize to find the best parameters
    println("Starting optimization...")

    result = bboptimize(
        cost; # cost functionw_
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)], # search space
        Method       = :adaptive_de_rand_1_bin_radiuslimited , #dxnes, or :cmaes, :genetic, etc.
        MaxSteps     = 1_000_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )
    println("Best solution = ", best_candidate(result))
    # 9) Extract the parameters
        ## Pumppulse
            ### dmin
                pulse_params_dmin  =  best_candidate(result)[1:3]
            ### rfr
                pulse_params_rfr  =  best_candidate(result)[4:6]
            ### rmin
                pulse_params_rmin  =  best_candidate(result)[7:9]
        ## Relaxation
            relaxation_params =  best_candidate(result)[10:13]
        ## Couplings
            kopplungs_params = best_candidate(result)[14:26]
        ## Degeneracy
            g = [1,2,12,16,1,3,2] # hardcoded for this model_rmin
    # 10) Create the the fitted model

    best_model_dmin  = model3_3(
        t_model_dmin,
        vcat(pulse_params_dmin,relaxation_params,kopplungs_params,g),
        mode=:dmin,
        dt=dt
    )
    best_model_rfr  = model3_3(
        t_model_rfr,
        vcat(pulse_params_rfr,relaxation_params,kopplungs_params,g),
        mode=:rfr,
        dt=dt
    )
    best_model_rmin  = model3_3(
        t_model_rmin,
        vcat(pulse_params_rmin,relaxation_params,kopplungs_params,g),
        mode=:rmin,
        dt=dt
    )

    return best_model_dmin,best_model_rfr,best_model_rmin
end


"""
    Modell 3.2
    From wavelength-scan experiments, it became clear that pumping at the wavelength labeled "rmin" actually excites two modes simultaneously 
    (rmina and rminb). Therefore, an additional pump pulse parameter has been introduced. This parameter is set to zero whenever "rmin" is not 
    actively pumped.
    Model 3.2 is considered an extension and refinement of Model 3.1.
    Modell 3.1:
    Model 3.1 is derived from Model 2.3. The primary shortcoming of Model 2.3 was its inability to accurately reproduce the fast dynamics 
    observed in the "rfr" mode. Model 3.1 addresses this issue by allowing forward and backward rates involving the "rfr" mode to differ (asymmetric transitions).
    As in Model 2.3, parameters p[1:2] describe the laser excitation processes. Parameters p[3:9] represent the relaxation rates of the modes, 
    whereas p[10:50] specify the coupling rates between different modes. Finally, parameters p[51:57] define the degeneracy of each mode.
        julia> model3_2(
                t,                  # model time vector
                p::Vector{T};       # model parameters
                mode = :dmin,       # which mode is pumped? :dmin, :rmin, or :rfr
                dt = 0.1            # time step size for the simulation
               ) where T<:Number

"""
function model3_2(
    t,
    p::Vector{T};
    mode = :dmin,
    dt = 0.1
) where T<:Number

    #---------------------------------------------------------
    # 1) We have 7 excited states. We'll store their populations
    #    in Popu[:,1] through Popu[:,7].
    #    We'll also store a "bleach" signal in Bleach[:,1..7].
    #---------------------------------------------------------
    n = 7  # number of states
    Popu   = zeros(T, length(t), n)
    Bleach = zeros(T, length(t), n)
    control = zeros(T, length(t), 4)

    #---------------------------------------------------------
    # 2) Unpack parameters from p:
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    α_rb     = 0      # amplitude scaling of laser pulse for rminb/u
    if mode == :rmin
        α_rb = p[2]
    end
    t_0      = p[3]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width to observed pulse width of pump pulse

    # Relaxation rates to ground state 
    relax_rplus       = 1/p[4]
    relax_dminusomega = 1/p[5]
    relax_dfr         = 1/p[6]
    relax_dminus      = 1/p[7]
    relax_rfr         = 1/p[8]
    relax_rminusb     = 1/p[9]
    relax_rminusa     = 1/p[10]

    # Coupling constants between states
    # For simplicity, I'll rename them 'k_x_y' to emphasize
    # they are coupling rates. Each will get multiplied by
    # the product of degeneracies g_x*g_y for symmetrical flow.

    k_rplus_dminusomega   = 1/p[11]
    k_rplus_dfr           = 1/p[12]
    k_rplus_dminus        = 1/p[13]
    k_rplus_rfr           = 1/p[14]
    k_rplus_rminusb       = 1/p[15]
    k_rplus_rminusa       = 1/p[16]

    k_dminusomega_rplus   = k_rplus_dminusomega
    k_dminusomega_dfr     = 1/p[17]
    k_dminusomega_dminus  = 1/p[18]
    k_dminusomega_rfr     = 1/p[19]
    k_dminusomega_rminusb = 1/p[20]
    k_dminusomega_rminusa = 1/p[21]

    k_dfr_rplus           = k_rplus_dfr
    k_dfr_dminusomega     = k_dminusomega_dfr
    k_dfr_dminus          = 1/p[22]
    k_dfr_rfr             = 1/p[23]
    k_dfr_rminusb         = 1/p[24]
    k_dfr_rminusa         = 1/p[25]

    k_dminus_rplus        = k_rplus_dminus
    k_dminus_dminusomega  = k_dminusomega_dminus
    k_dminus_dfr          = k_dfr_dminus
    k_dminus_rfr          = 1/p[26]
    k_dminus_rminusb      = 1/p[27]
    k_dminus_rminusa      = 1/p[28]

    k_rfr_rplus           = 1/p[29]
    k_rfr_dminusomega     = 1/p[30]
    k_rfr_dfr             = 1/p[31]
    k_rfr_dminus          = 1/p[32]
    k_rfr_rminusb         = 1/p[33]
    k_rfr_rminusa         = 1/p[34]

    k_rminusb_rplus       = k_rplus_rminusb
    k_rminusb_dminusomega = k_dminusomega_rminusb
    k_rminusb_dfr         = k_dfr_rminusb
    k_rminusb_dminus      = k_dminus_rminusb
    k_rminusb_rfr         = 1/p[35]
    k_rminusb_rminusa     = 1/p[36]

    k_rminusa_rplus       = k_rplus_rminusa
    k_rminusa_dminusomega = k_dminusomega_rminusa
    k_rminusa_dfr         = k_dfr_rminusa
    k_rminusa_dminus      = k_dminus_rminusa
    k_rminusa_rfr         = 1/p[37]
    k_rminusa_rminusb     = k_rminusb_rminusa

    # 3) Degeneracy of the states:
    # only fixed values for rplus and rfr
    g_rplus   = 1 #p[38]
    g_d_omega = round(Int,p[39])
    g_dfr     = round(Int,p[40])
    g_dminus  = round(Int,p[41])
    g_rfr     = 1 #p[42]
    g_rminb   = round(Int,p[43])
    g_rmina   = round(Int,p[44])


    # Store them in an array:
    #   1->rplus, 2->dminusomega, 3->dfr, 4->dminus, 
    #   5->rfr,   6->rminusb,     7->rminusa
        G = [g_rplus, g_d_omega, g_dfr, g_dminus, g_rfr, g_rminb, g_rmina]

    #---------------------------------------------------------
    # 4) Precompute the pulse shape at each time point
    f_g    = gauss.(t,Ref([σ,t_0]))

    #---------------------------------------------------------
    # 5) Mode-specific initial conditions:
    #    :dmin, :rmin, or :rfr  -> which state is initially pumped?
    #---------------------------------------------------------
        pump_first = α * f_g[1] * dt
    # Initially, one state is pumped:
        Popu[1,1] = (mode == :rplus)     ? 1/G[1] * pump_first         : 0
        Popu[1,2] = (mode == :dminomega) ? 1/G[2] * pump_first         : 0
        Popu[1,3] = (mode == :dfr)       ? 1/G[3] * pump_first         : 0
        Popu[1,4] = (mode == :dmin)      ? 1/G[4] * pump_first         : 0
        Popu[1,5] = (mode == :rfr)       ? 1/G[5] * pump_first         : 0
        Popu[1,6] = (mode == :rmin)      ? 1/G[6] * α_rb * f_g[1] * dt : 0
        Popu[1,7] = (mode == :rmin)      ? 1/G[7] * pump_first         : 0
        control[1,1] = pump_first + α_rb * f_g[1] * dt 
        control[1,2] = G[1] * Popu[1,1] + (
            + G[2] * Popu[1,2] 
            + G[3] * Popu[1,3] 
            + G[4] * Popu[1,4] 
            + G[5] * Popu[1,5] 
            + G[6] * Popu[1,6] 
            + G[7] * Popu[1,7]
        )


    #---------------------------------------------------------
    # 6) Time stepping with explicit Euler method:
    #    We'll update all 7 states at each time step.
    #
    #    Key change for symmetrical coupling:
    #      * For a coupling k_x_y between state x and y,
    #        we multiply by (G[x]*G[y]) in BOTH equations:
    #
    #        Popu[i,x] += + k_xy * G[x]*G[y]*(Popu[i-1,y] - Popu[i-1,x])
    #        Popu[i,y] += + k_xy * G[x]*G[y]*(Popu[i-1,x] - Popu[i-1,y])
    #
    #    The same factor appears in both directions.
    #---------------------------------------------------------
    for i in 2:length(t)

        # Precompute the laser amplitude at this timestep
        laser_here = α * f_g[i]

        #-------------
        # State 1: rplus 
        #-------------
        # If this is the pumped state in :rplus mode, we add laser_here
        # (same code runs for all modes, but laser_here=0 if not used)
        pump_rplus = (mode == :rplus) ?  laser_here : 0
        Popu[i,1] = Popu[i-1,1] + 1/G[1] * (

            # Pumping if this mode is :rplus
            + pump_rplus

            # Relaxation to ground
            - G[1] * relax_rplus * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rplus_dminusomega * (G[1]*G[2]) * Popu[i-1,1]
            + k_dminusomega_rplus * (G[2]*G[1]) * Popu[i-1,2] 

            # Coupling with dfr (state 3)
            - k_rplus_dfr         * (G[1]*G[3]) * Popu[i-1,1]
            + k_dfr_rplus         * (G[3]*G[1]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rplus_dminus      * (G[1]*G[4]) * Popu[i-1,1]
            + k_dminus_rplus      * (G[4]*G[1]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rplus_rfr         * (G[1]*G[5]) * Popu[i-1,1]
            + k_rfr_rplus         * (G[5]*G[1]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rplus_rminusb     * (G[1]*G[6]) * Popu[i-1,1]
            + k_rminusb_rplus     * (G[6]*G[1]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rplus_rminusa     * (G[1]*G[7]) * Popu[i-1,1]
            + k_rminusa_rplus     * (G[7]*G[1]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 2: dminusomega
        #-------------
        pump_dminomega = (mode == :dminomega) ? laser_here : 0
        Popu[i,2] = Popu[i-1,2] + 1/G[2] * (

            # Pumping if this mode is :dminomega
            + pump_dminomega

            # Relaxation
            - G[2] * relax_dminusomega * Popu[i-1,2]

            # Coupling with rplus (state 1)
            - k_dminusomega_rplus   * (G[2]*G[1]) * Popu[i-1,2]
            + k_rplus_dminusomega   * (G[1]*G[2]) * Popu[i-1,1]

            # Coupling with dfr (state 3)
            - k_dminusomega_dfr     * (G[2]*G[3]) * Popu[i-1,2]
            + k_dfr_dminusomega     * (G[3]*G[2]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_dminusomega_dminus  * (G[2]*G[4]) * Popu[i-1,2]
            + k_dminus_dminusomega  * (G[4]*G[2]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dminusomega_rfr     * (G[2]*G[5]) * Popu[i-1,2]
            + k_rfr_dminusomega     * (G[5]*G[2]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]
            + k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]
            + k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 3: dfr
        #-------------
        pump_dfr = (mode == :dfr) ? laser_here : 0
        Popu[i,3] = Popu[i-1,3] + 1/G[3] * (

            # Pumping if this mode is :dfr
            + pump_dfr

            # Relaxation
            - G[3] * relax_dfr * Popu[i-1,3]

            # Coupling with rplus (state 1)
            - k_dfr_rplus       * (G[3]*G[1]) * Popu[i-1,3]
            + k_rplus_dfr       * (G[1]*G[3]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dfr_dminusomega * (G[3]*G[2]) * Popu[i-1,3]
            + k_dminusomega_dfr * (G[2]*G[3]) * Popu[i-1,2]

            # Coupling with dminus (state 4)
            - k_dfr_dminus      * (G[3]*G[4]) * Popu[i-1,3]
            + k_dminus_dfr      * (G[4]*G[3]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]
            + k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dfr_rminusb     * (G[3]*G[6]) * Popu[i-1,3]
            + k_rminusb_dfr     * (G[6]*G[3]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dfr_rminusa     * (G[3]*G[7]) * Popu[i-1,3]
            + k_rminusa_dfr     * (G[7]*G[3]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 4: dminus
        #-------------
        pump_dmin = (mode == :dmin) ? laser_here : 0
        Popu[i,4] = Popu[i-1,4] + 1/G[4] * (

            # Pumping if this mode is :dmin
            + pump_dmin

            # Relaxation
            - G[4] * relax_dminus * Popu[i-1,4]

            # Coupling with rplus (state 1)
            - k_dminus_rplus       * (G[4]*G[1]) * Popu[i-1,4]
            + k_rplus_dminus       * (G[1]*G[4]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dminus_dminusomega * (G[4]*G[2]) * Popu[i-1,4]
            + k_dminusomega_dminus * (G[2]*G[4]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_dminus_dfr         * (G[4]*G[3]) * Popu[i-1,4]
            + k_dfr_dminus         * (G[3]*G[4]) * Popu[i-1,3]

            # Coupling with rfr (state 5)
            - k_dminus_rfr         * (G[4]*G[5]) * Popu[i-1,4]
            + k_rfr_dminus         * (G[5]*G[4]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminus_rminusb     * (G[4]*G[6]) * Popu[i-1,4]
            + k_rminusb_dminus     * (G[6]*G[4]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminus_rminusa     * (G[4]*G[7]) * Popu[i-1,4]
            + k_rminusa_dminus     * (G[7]*G[4]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 5: rfr
        #-------------
        pump_rfr = (mode == :rfr) ? laser_here : 0
        Popu[i,5] = Popu[i-1,5] + 1/G[5] * (
            # Pumping if this mode is :rfr
            + pump_rfr
            
            # Relaxation
            - G[5] * relax_rfr * Popu[i-1,5]

            # Coupling with rplus (state 1)
            - k_rfr_rplus       * (G[5]*G[1]) * Popu[i-1,5]
            + k_rplus_rfr       * (G[1]*G[5]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rfr_dminusomega * (G[5]*G[2]) * Popu[i-1,5]
            + k_dminusomega_rfr * (G[2]*G[5]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]
            + k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rfr_dminus      * (G[5]*G[4]) * Popu[i-1,5]
            + k_dminus_rfr      * (G[4]*G[5]) * Popu[i-1,4]

            # Coupling with rminusb (state 6)
            - k_rfr_rminusb     * (G[5]*G[6]) * Popu[i-1,5]
            + k_rminusb_rfr     * (G[6]*G[5]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rfr_rminusa     * (G[5]*G[7]) * Popu[i-1,5]
            + k_rminusa_rfr     * (G[7]*G[5]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 6: rminusb
        #-------------
        #pump_rminb = (mode == :rminb) ? laser_here : 0
        pump_rminb  = α_rb * f_g[i] * dt
        Popu[i,6] = Popu[i-1,6] + 1/G[6] * (

            # Pumping if this mode is :rminb
            + pump_rminb

            # Relaxation
            - G[6] * relax_rminusb * Popu[i-1,6]

            # Coupling with rplus (state 1)
            - k_rminusb_rplus       * (G[6]*G[1]) * Popu[i-1,6]
            + k_rplus_rminusb       * (G[1]*G[6]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]
            + k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusb_dfr         * (G[6]*G[3]) * Popu[i-1,6]
            + k_dfr_rminusb         * (G[3]*G[6]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusb_dminus      * (G[6]*G[4]) * Popu[i-1,6]
            + k_dminus_rminusb      * (G[4]*G[6]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusb_rfr         * (G[6]*G[5]) * Popu[i-1,6]
            + k_rfr_rminusb         * (G[5]*G[6]) * Popu[i-1,5]

            # Coupling with rminusa (state 7)
            - k_rminusb_rminusa     * (G[6]*G[7]) * Popu[i-1,6]
            + k_rminusa_rminusb     * (G[7]*G[6]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 7: rminusa
        #-------------
        pump_rmin = (mode == :rmin) ? laser_here : 0
        Popu[i,7] = Popu[i-1,7] + 1/G[7] * (

            # Pumping if this mode is :rmin
            + pump_rmin

            # Relaxation
            - G[7] * relax_rminusa * Popu[i-1,7]

            # Coupling with rplus (state 1)
            - k_rminusa_rplus       * (G[7]*G[1]) * Popu[i-1,7]
            + k_rplus_rminusa       * (G[1]*G[7]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
            + k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusa_dfr         * (G[7]*G[3]) * Popu[i-1,7]
            + k_dfr_rminusa         * (G[3]*G[7]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusa_dminus      * (G[7]*G[4]) * Popu[i-1,7]
            + k_dminus_rminusa      * (G[4]*G[7]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusa_rfr        * (G[7]*G[5]) * Popu[i-1,7]
            + k_rfr_rminusa        * (G[5]*G[7]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rminusa_rminusb    * (G[7]*G[6]) * Popu[i-1,7]
            + k_rminusb_rminusa    * (G[6]*G[7]) * Popu[i-1,6]
        ) * dt

        #-------------
        # 7) Control
        # Checks if the balance of the populations is correct at each timestep
        #-------------

        # 1. Integral of Laser Input 
        control[i,1] = control[i-1,1] + (laser_here + α_rb * f_g[i]) * dt
        # 2. Overall Population of all substates 
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )
        
        # 3. Relaxed Population of all substates
        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        
        # 4. Sum of Laser Input, Overall Population and Relaxed Population should equal to 0 at every timestep
        control[i,4] = control[i,1] - control[i,2] - control[i,3]
    end

    #---------------------------------------------------------
    # 8) Finally, compute the "Bleach" signal 
    #    for each state: (1 - 2*g[i]*Popu[i])^2
    #---------------------------------------------------------
    for k in 1:n
        Bleach[:,k] = (1 .- 2 .* G[k] .* Popu[:,k]) .^2
    end

# ------------------------------------------------------
    # 9) Dictionary (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1],alpha_rb = p[2], t_0 = p[3])

    # B) Relaxation
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[4:10]
    )

    # C) Entartung
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G#p[38:44]
    )

    # D) Couplings between states
    #    We differ between forward and backward couplings
    #    Syntax: (x, y, fwd_idx, bwd_idx)
    #      => t_x_y   => p[fwd_idx]
    #      => t_y_x   => p[bwd_idx]
    couplings = [
        (:rplus,       :dminusomega, 11, 11),
        (:rplus,       :dfr,         12, 12),
        (:rplus,       :dminus,      13, 13),
        (:rplus,       :rfr,         14, 29),  # asym
        (:rplus,       :rminusb,     15, 15),
        (:rplus,       :rminusa,     16, 16),

        (:dminusomega, :dfr,         17, 17),
        (:dminusomega, :dminus,      18, 18),
        (:dminusomega, :rfr,         19, 30),  # asym
        (:dminusomega, :rminusb,     20, 20),
        (:dminusomega, :rminusa,     21, 21),

        (:dfr,         :dminus,      22, 22),
        (:dfr,         :rfr,         23, 31),  # asym
        (:dfr,         :rminusb,     24, 24),
        (:dfr,         :rminusa,     25, 25),

        (:dminus,      :rfr,         26, 32),  # asym
        (:dminus,      :rminusb,     27, 27),
        (:dminus,      :rminusa,     28, 28),

        (:rfr,         :rminusb,     33, 35),  # asym
        (:rfr,         :rminusa,     34, 37),  # asym

        (:rminusb,     :rminusa,     36, 36)
    ]

    rows_all = Vector{Dict{Symbol,Any}}()
    for (x, y, fwd_idx, bwd_idx) in couplings
        push!(rows_all, Dict(
            :param_name => "t_$(x)_$(y)",
            :value      => p[fwd_idx],
            :ratio     => p[fwd_idx]/p[bwd_idx]

        ))
        push!(rows_all, Dict(
            :param_name => "t_$(y)_$(x)",
            :value      => p[bwd_idx],
            :ratio     => p[bwd_idx]/p[fwd_idx]
        ))
    end
    df_koppl_all = DataFrame(rows_all)

    # Function to create DataFrame for each mode
    function df_koppl_for_mode(m::Symbol)
        subrows = Vector{Dict{Symbol,Any}}()
        for (x,y, fwd_idx, bwd_idx) in couplings
            if x == m || y == m
                # x->y
                push!(subrows, Dict(
                    :param_name => "t_$(x)_$(y)",
                    :value      => p[fwd_idx],
                    :ratio     => p[fwd_idx]/p[bwd_idx]
                ))
                # y->x
                push!(subrows, Dict(
                    :param_name => "t_$(y)_$(x)",
                    :value      => p[bwd_idx],
                    :ratio     => p[bwd_idx]/p[fwd_idx]
                ))
            end
        end
        return DataFrame(subrows)
    end

    # Create DataFrame for each mode
    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Full output dictionary
    outputDict = Dict(
        :Model       => "Model 3.2",
        :fitted_mode => mode,
        :Parameter   => Dict(
            :All        => p,
            :Puls       => df_puls,
            :Relaxation => df_relax,
            :Entartung  => df_entartung,
            :Kopplung   => dict_koppl
        ),
        :Zeitschritt => dt
    )

    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end

"""
    The following function can fit the model 3.2 to the data. For more information on the model, see the model function model3_2.
    The function uses a global optimizer of the package "BlackBoxOptim.jl". It is recommend to use the :adaptive_de_rand_1_bin_radiuslimited algorithm 
    to have a robustfree fit. 
    This function fits three data sets to one set of parameters. The data sets are:
        - time_dmin, traces_dmin, use_dmin: data set for the dmin pumped experiment
        - time_rfr, traces_rfr, use_rfr: data set for the rfr pumped experiment
        - time_rmin, traces_rmin, use_rmin: data set for the rmin pumped experiment
    The time argument is the time vector of the corresponding data set, the traces argument is the data set itself 
    of the used traces and the use argument is a vector which contains the indices of the traces which should be used for the fit.
    The coha argument is an array. coha[1] contains the index of the mode where we observed the coherent artifact in the rfr pumped experiment.
    coha[2] contains the indices of the time_vector where the coherent artifact is observed.
    The function returns three ReservoirModel objects for the three data sets.
        julia> blackbox_model_3_2(
                time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
                time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
                time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
                coha;                           # cohrent artifact array
                dt=0.1                          # time step
               )
"""
function blackbox_model_3_2(
    time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
    time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
    time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
    coha;                           # cohrent artifact array
    dt=0.1                          # time step
)
    # 1) Define the search space
    lb = vcat(
        [0.0, 0.0, -20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [0.0, 0.0, -20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [0.0, 0.0, -20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(1,7), # relaxation_params   
        fill(1,27), # kopplungs_params
        1.00,  # g_rplus
        1.0,   # g_dminomega
        1.0,   # g_dfr
        8.0,   # g_dminus
        1.00,  # g_rfr
        1.00,  # g_rminusb
        1.99,  # g_rminusa 
    )
    ub = vcat(
        [50.0, 1e-5, 20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [50.0, 1e-5, 20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [50.0, 50.0, 20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(20000.0,7), # relaxation_params   
        fill(20000.0,27), # kopplungs_params
        1.01,              # g_rplus
        2.00,              # g_dminomega
        10.0,             # g_dfr
        10.0,             # g_dminus
        1.01,              # g_rfr
        3.0,              # g_rminusb
        2.01               # g_rminusa
    )
    # 2) Construct time vectors for each experiment
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    # 3) Define the cost function
    # The function "cost" computes a total cost based on the weighted squared differences
    # between the model predictions and experimental traces. It includes penalties to
    # constrain the simulated populations.
    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1]; α_rb_dmin  = p[2];  t0_dmin  = p[3]
        α_rfr   = p[4]; α_rb_rfr   = p[5];  t0_rfr   = p[6]
        α_rmin  = p[7]; α_rb_rmin  = p[8];  t0_rmin  = p[9]

        relax_p = p[10:16]      # 7 relaxation
        koppl_p = p[17:43]     # 27 couplings
        g     = p[44:end]      # 7 degeneracy

        # 2) Construct the model
        # With the use argument we select a subset of the model to fit to the data

        # dmin pumped experiment
        model_dmin = model3_2(
            t_model_dmin,
            vcat(α_dmin,α_rb_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # dmin subset
        model_dmin_used = model_dmin.Bleach[:, use_dmin]

        # rfr pumped experiment
        model_rfr = model3_2(
            t_model_rfr,
            vcat(α_rfr, α_rb_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # rfr subset
        model_rfr_used = model_rfr.Bleach[:, use_rfr]

        # rmin pumped experiment
        model_rmin = model3_2(
            t_model_rmin,
            vcat(α_rmin,α_rb_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # rmin subset
        model_rmin_used = model_rmin.Bleach[:, use_rmin]

        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
        #   Also, replace NaN and Inf values with a large number (1e12)
        pred_dmin = hcat(v_pred_dmin...)
        for i in 1:length(pred_dmin)
            if isinf(pred_dmin[i]) || isnan(pred_dmin[i])
                pred_dmin[i] = 1e12
            end
        end
        pred_rfr  = hcat(v_pred_rfr...)
        for i in 1:length(pred_rfr)
            if isinf(pred_rfr[i]) || isnan(pred_rfr[i])
                pred_rfr[i] = 1e12
            end
        end
        pred_rmin = hcat(v_pred_rmin...)
        for i in 1:length(pred_rmin)
            if isinf(pred_rmin[i]) || isnan(pred_rmin[i])
                pred_rmin[i] = 1e12
            end
        end
        # 5) Select time steps without coherence artifact in rfr pumped experiment
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights 
        # Use the time_weight function to create special weights for short delays, since we observe 
        # fast dynamics in the beginning of the experiment.

        w_dmin = ones(size(pred_dmin))
        w_pertrace_dmin = [
            5.0,  # rplus
            1.0,  # dfr
            5.0,  # rfr
            5.0   # rminusa

        ]
        for m in 1:size(w_dmin,1)
            for n in 1:size(w_dmin,2)
                w_dmin[m,n] = time_weight(time_dmin[m],weight=3) * w_pertrace_dmin[n]
            end
        end
    
        w_rfr = ones(size(pred_rfr))
        w_pertrace_rfr = [
            5.0,  # rplus
            1.0,  # dfr
            5.0,  # rfr
            5.0   # rminusa
        ]
        for m in 1:size(w_rfr,1)
            for n in 1:size(w_rfr,2)
                if (n == coha[1]) && (m ∉ selected_idx_rfr)
                    w_rfr[m,n] = 0.0
                else
                    w_rfr[m,n] = time_weight(time_rfr[m],weight=3) * w_pertrace_rfr[n]
                end
            end
        end
        w_rmin = ones(size(pred_rmin))
        w_pertrace_rmin = [
            5.0,  # rplus
            5.0,  # rfr
            2.0,  # rminusb
            5.0   # rminusa
        ]
        for m in 1:size(w_rmin,1)
            for n in 1:size(w_rmin,2)
                if n == 9
                    w_rmin[m,n] = 10 .* w_pertrace_rmin[n]
                else
                    w_rmin[m,n] = time_weight(time_rmin[m],weight=3) * w_pertrace_rmin[n]
                end
            end
        end

        # 7) Since we sum over all Data points, we need to normalize the weights
        # to the length of the given dataset
            length_dmin = length(time_dmin)
            length_rfr = length(time_rfr)
            length_rmin = length(time_rmin)
            norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
            norm_weight ./= maximum(norm_weight)
            cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
            cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
            cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 8) Penaltys
        # We define a maximum Population with max_P
        # If the Population is above this value, we add a penalty to the cost function
            max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
            P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
                norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
                norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
            total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty


        return total_cost
    end

    # Use bboptimize to find the best parameters
    println("Starting optimization...")

    result = bboptimize(
        cost; # cost function
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)], # search space
        Method       = :adaptive_de_rand_1_bin_radiuslimited ,   # or :cmaes, :genetic, etc.
        MaxSteps     = 10_000_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )

    println("Best solution = ", best_candidate(result))
    # 9) Extract the parameters
        ## Pumppulse
            ### dmin
                pulse_params_dmin  =  best_candidate(result)[1:3]
            ### rfr
                pulse_params_rfr  =  best_candidate(result)[4:6]
            ### rmin
                pulse_params_rmin  =  best_candidate(result)[7:9]
        ## Relaxation
            relaxation_params =  best_candidate(result)[10:16]
        ## Kopplung
            kopplungs_params = best_candidate(result)[17:43]
        ## Degeneracy
            g = best_candidate(result)[44:end]
    # 10) Create the the fitted model

    best_model_dmin  = model3_2(
        t_model_dmin,
        vcat(pulse_params_dmin,relaxation_params,kopplungs_params,g),
        mode=:dmin,
        dt=dt
    )
    best_model_rfr  = model3_2(
        t_model_rfr,
        vcat(pulse_params_rfr,relaxation_params,kopplungs_params,g),
        mode=:rfr,
        dt=dt
    )
    best_model_rmin  = model3_2(
        t_model_rmin,
        vcat(pulse_params_rmin,relaxation_params,kopplungs_params,g),
        mode=:rmin,
        dt=dt
    )

    return best_model_dmin,best_model_rfr,best_model_rmin
end

"""
    Modell 3.1:
    Model 3.1 is derived from Model 2.3. The primary shortcoming of Model 2.3 was its inability to accurately reproduce the fast dynamics 
    observed in the "rfr" mode. Model 3.1 addresses this issue by allowing forward and backward rates involving the "rfr" mode to differ (asymmetric transitions).
    As in Model 2.3, parameters p[1:2] describe the laser excitation processes. Parameters p[3:9] represent the relaxation rates of the modes, 
    whereas p[10:50] specify the coupling rates between different modes. Finally, parameters p[51:57] define the degeneracy of each mode.
        julia> model3_1(
                t,                  # model time vector
                p::Vector{T};       # model parameters
                mode = :dmin,       # which mode is pumped? :dmin, :rmin, or :rfr
                dt = 0.1            # time step size for the simulation
               ) where T<:Number

"""
function model3_1(
    t,
    p::Vector{T};
    mode = :dmin,
    dt = 0.1
) where T<:Number

    #---------------------------------------------------------
    # 1) We have 7 excited states. We'll store their populations
    #    in Popu[:,1] through Popu[:,7].
    #    We'll also store a "bleach" signal in Bleach[:,1..7].
    #---------------------------------------------------------
    n = 7  # number of states
    Popu   = zeros(T, length(t), n)
    Bleach = zeros(T, length(t), n)
    control = zeros(T, length(t), 4)

    #---------------------------------------------------------
    # 2) Unpack parameters from p:
    #    p[1], p[2], etc.  (Keep same indexing as your code.)
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    t_0      = p[2]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width to observed pulse width of pump pulse

    # Relaxation rates to ground state 
    relax_rplus       = 1/p[3]
    relax_dminusomega = 1/p[4]
    relax_dfr         = 1/p[5]
    relax_dminus      = 1/p[6]
    relax_rfr         = 1/p[7]
    relax_rminusb     = 1/p[8]
    relax_rminusa     = 1/p[9]

    # Coupling constants between states 
    # For simplicity, I'll rename them 'k_x_y' to emphasize
    # they are coupling rates. Each will get multiplied by
    # the product of degeneracies g_x*g_y for symmetrical flow.

    k_rplus_dminusomega   = 1/p[10]
    k_rplus_dfr           = 1/p[11]
    k_rplus_dminus        = 1/p[12]
    k_rplus_rfr           = 1/p[13]
    k_rplus_rminusb       = 1/p[14]
    k_rplus_rminusa       = 1/p[15]

    k_dminusomega_rplus   = k_rplus_dminusomega
    k_dminusomega_dfr     = 1/p[16]
    k_dminusomega_dminus  = 1/p[17]
    k_dminusomega_rfr     = 1/p[18]
    k_dminusomega_rminusb = 1/p[19]
    k_dminusomega_rminusa = 1/p[20]

    k_dfr_rplus           = k_rplus_dfr
    k_dfr_dminusomega     = k_dminusomega_dfr
    k_dfr_dminus          = 1/p[21]
    k_dfr_rfr             = 1/p[22]
    k_dfr_rminusb         = 1/p[23]
    k_dfr_rminusa         = 1/p[24]

    k_dminus_rplus        = k_rplus_dminus
    k_dminus_dminusomega  = k_dminusomega_dminus
    k_dminus_dfr          = k_dfr_dminus
    k_dminus_rfr          = 1/p[25]
    k_dminus_rminusb      = 1/p[26]
    k_dminus_rminusa      = 1/p[27]

    k_rfr_rplus           = 1/p[28]
    k_rfr_dminusomega     = 1/p[29]
    k_rfr_dfr             = 1/p[30]
    k_rfr_dminus          = 1/p[31]
    k_rfr_rminusb         = 1/p[32]
    k_rfr_rminusa         = 1/p[33]

    k_rminusb_rplus       = k_rplus_rminusb
    k_rminusb_dminusomega = k_dminusomega_rminusb
    k_rminusb_dfr         = k_dfr_rminusb
    k_rminusb_dminus      = k_dminus_rminusb
    k_rminusb_rfr         = 1/p[34]
    k_rminusb_rminusa     = 1/p[35]

    k_rminusa_rplus       = k_rplus_rminusa
    k_rminusa_dminusomega = k_dminusomega_rminusa
    k_rminusa_dfr         = k_dfr_rminusa
    k_rminusa_dminus      = k_dminus_rminusa
    k_rminusa_rfr         = 1/p[36]
    k_rminusa_rminusb     = k_rminusb_rminusa

    # 3) Degeneracy of the states:
    # Fixed values after Model 3.2
    g_rplus   = p[37]
    g_d_omega = p[38]
    g_dfr     = p[39]
    g_dminus  = p[40]
    g_rfr     = p[41]
    g_rminb   = p[42]
    g_rmina   = p[43]


    # Store them in an array:
    #   1->rplus, 2->dminusomega, 3->dfr, 4->dminus, 
    #   5->rfr,   6->rminusb,     7->rminusa
    G = [g_rplus, g_d_omega, g_dfr, g_dminus, g_rfr, g_rminb, g_rmina]

    #---------------------------------------------------------
    # 4) Precompute the pulse shape at each time point
    f_g    = gauss.(t,Ref([σ,t_0]))

    #---------------------------------------------------------
    # 5) Mode-specific initial conditions:
    #    :dmin, :rmin, or :rfr  -> which state is initially pumped?
    #---------------------------------------------------------
        pump_first = α * f_g[1] * dt
    # Initially, one state is pumped:
        Popu[1,1] = (mode == :rplus)     ? 1/G[1] * pump_first : 0
        Popu[1,2] = (mode == :dminomega) ? 1/G[2] * pump_first : 0
        Popu[1,3] = (mode == :dfr)       ? 1/G[3] * pump_first : 0
        Popu[1,4] = (mode == :dmin)      ? 1/G[4] * pump_first : 0
        Popu[1,5] = (mode == :rfr)       ? 1/G[5] * pump_first : 0
        Popu[1,6] = (mode == :rminb)     ? 1/G[6] * pump_first : 0
        Popu[1,7] = (mode == :rmin)      ? 1/G[7] * pump_first : 0
        control[1,1] = pump_first
        control[1,2] = G[1] * Popu[1,1] + (
            + G[2] * Popu[1,2] 
            + G[3] * Popu[1,3] 
            + G[4] * Popu[1,4] 
            + G[5] * Popu[1,5] 
            + G[6] * Popu[1,6] 
            + G[7] * Popu[1,7]
        )


    #---------------------------------------------------------
    # 6) Time stepping with explicit Euler method:
    #    We'll update all 7 states at each time step.
    #
    #    Key change for symmetrical coupling:
    #      * For a coupling k_xy between state x and y,
    #        we multiply by (G[x]*G[y]) in BOTH equations:
    #
    #        Popu[i,x] += + k_xy * G[x]*G[y]*(Popu[i-1,y] - Popu[i-1,x])
    #        Popu[i,y] += + k_xy * G[x]*G[y]*(Popu[i-1,x] - Popu[i-1,y])
    #
    #    The same factor appears in both directions.
    #---------------------------------------------------------
    for i in 2:length(t)

        # Precompute the laser amplitude at this time, if needed
        laser_here = α * f_g[i]

        #-------------
        # State 1: rplus 
        #-------------
        # If this is the pumped state in :rplus mode, we add laser_here
        # (same code runs for all modes, but laser_here=0 if not used)
        pump_rplus = (mode == :rplus) ?  laser_here : 0
        Popu[i,1] = Popu[i-1,1] + 1/G[1] * (

            # Pumping if this mode is :rplus
            + pump_rplus

            # Relaxation to ground
            - G[1] * relax_rplus * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rplus_dminusomega * (G[1]*G[2]) * Popu[i-1,1]
            + k_dminusomega_rplus * (G[2]*G[1]) * Popu[i-1,2] 

            # Coupling with dfr (state 3)
            - k_rplus_dfr         * (G[1]*G[3]) * Popu[i-1,1]
            + k_dfr_rplus         * (G[3]*G[1]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rplus_dminus      * (G[1]*G[4]) * Popu[i-1,1]
            + k_dminus_rplus      * (G[4]*G[1]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rplus_rfr         * (G[1]*G[5]) * Popu[i-1,1]
            + k_rfr_rplus         * (G[5]*G[1]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rplus_rminusb     * (G[1]*G[6]) * Popu[i-1,1]
            + k_rminusb_rplus     * (G[6]*G[1]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rplus_rminusa     * (G[1]*G[7]) * Popu[i-1,1]
            + k_rminusa_rplus     * (G[7]*G[1]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 2: dminusomega
        #-------------
        pump_dminomega = (mode == :dminomega) ? laser_here : 0
        Popu[i,2] = Popu[i-1,2] + 1/G[2] * (

            # Pumping if this mode is :dminomega
            + pump_dminomega

            # Relaxation
            - G[2] * relax_dminusomega * Popu[i-1,2]

            # Coupling with rplus (state 1)
            - k_dminusomega_rplus   * (G[2]*G[1]) * Popu[i-1,2]
            + k_rplus_dminusomega   * (G[1]*G[2]) * Popu[i-1,1]

            # Coupling with dfr (state 3)
            - k_dminusomega_dfr     * (G[2]*G[3]) * Popu[i-1,2]
            + k_dfr_dminusomega     * (G[3]*G[2]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_dminusomega_dminus  * (G[2]*G[4]) * Popu[i-1,2]
            + k_dminus_dminusomega  * (G[4]*G[2]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dminusomega_rfr     * (G[2]*G[5]) * Popu[i-1,2]
            + k_rfr_dminusomega     * (G[5]*G[2]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]
            + k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]
            + k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 3: dfr
        #-------------
        pump_dfr = (mode == :dfr) ? laser_here : 0
        Popu[i,3] = Popu[i-1,3] + 1/G[3] * (

            # Pumping if this mode is :dfr
            + pump_dfr

            # Relaxation
            - G[3] * relax_dfr * Popu[i-1,3]

            # Coupling with rplus (state 1)
            - k_dfr_rplus       * (G[3]*G[1]) * Popu[i-1,3]
            + k_rplus_dfr       * (G[1]*G[3]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dfr_dminusomega * (G[3]*G[2]) * Popu[i-1,3]
            + k_dminusomega_dfr * (G[2]*G[3]) * Popu[i-1,2]

            # Coupling with dminus (state 4)
            - k_dfr_dminus      * (G[3]*G[4]) * Popu[i-1,3]
            + k_dminus_dfr      * (G[4]*G[3]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]
            + k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dfr_rminusb     * (G[3]*G[6]) * Popu[i-1,3]
            + k_rminusb_dfr     * (G[6]*G[3]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dfr_rminusa     * (G[3]*G[7]) * Popu[i-1,3]
            + k_rminusa_dfr     * (G[7]*G[3]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 4: dminus
        #-------------
        pump_dmin = (mode == :dmin) ? laser_here : 0
        Popu[i,4] = Popu[i-1,4] + 1/G[4] * (

            # Pumping if this mode is :dmin
            + pump_dmin

            # Relaxation
            - G[4] * relax_dminus * Popu[i-1,4]

            # Coupling with rplus (state 1)
            - k_dminus_rplus       * (G[4]*G[1]) * Popu[i-1,4]
            + k_rplus_dminus       * (G[1]*G[4]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dminus_dminusomega * (G[4]*G[2]) * Popu[i-1,4]
            + k_dminusomega_dminus * (G[2]*G[4]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_dminus_dfr         * (G[4]*G[3]) * Popu[i-1,4]
            + k_dfr_dminus         * (G[3]*G[4]) * Popu[i-1,3]

            # Coupling with rfr (state 5)
            - k_dminus_rfr         * (G[4]*G[5]) * Popu[i-1,4]
            + k_rfr_dminus         * (G[5]*G[4]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminus_rminusb     * (G[4]*G[6]) * Popu[i-1,4]
            + k_rminusb_dminus     * (G[6]*G[4]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminus_rminusa     * (G[4]*G[7]) * Popu[i-1,4]
            + k_rminusa_dminus     * (G[7]*G[4]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 5: rfr
        #-------------
        pump_rfr = (mode == :rfr) ? laser_here : 0
        Popu[i,5] = Popu[i-1,5] + 1/G[5] * (
            # Pumping if this mode is :rfr
            + pump_rfr
            
            # Relaxation
            - G[5] * relax_rfr * Popu[i-1,5]

            # Coupling with rplus (state 1)
            - k_rfr_rplus       * (G[5]*G[1]) * Popu[i-1,5]
            + k_rplus_rfr       * (G[1]*G[5]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rfr_dminusomega * (G[5]*G[2]) * Popu[i-1,5]
            + k_dminusomega_rfr * (G[2]*G[5]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]
            + k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rfr_dminus      * (G[5]*G[4]) * Popu[i-1,5]
            + k_dminus_rfr      * (G[4]*G[5]) * Popu[i-1,4]

            # Coupling with rminusb (state 6)
            - k_rfr_rminusb     * (G[5]*G[6]) * Popu[i-1,5]
            + k_rminusb_rfr     * (G[6]*G[5]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rfr_rminusa     * (G[5]*G[7]) * Popu[i-1,5]
            + k_rminusa_rfr     * (G[7]*G[5]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 6: rminusb
        #-------------
        pump_rminb = (mode == :rminb) ? laser_here : 0
        Popu[i,6] = Popu[i-1,6] + 1/G[6] * (

            # Pumping if this mode is :rminb
            + pump_rminb

            # Relaxation
            - G[6] * relax_rminusb * Popu[i-1,6]

            # Coupling with rplus (state 1)
            - k_rminusb_rplus       * (G[6]*G[1]) * Popu[i-1,6]
            + k_rplus_rminusb       * (G[1]*G[6]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]
            + k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusb_dfr         * (G[6]*G[3]) * Popu[i-1,6]
            + k_dfr_rminusb         * (G[3]*G[6]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusb_dminus      * (G[6]*G[4]) * Popu[i-1,6]
            + k_dminus_rminusb      * (G[4]*G[6]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusb_rfr         * (G[6]*G[5]) * Popu[i-1,6]
            + k_rfr_rminusb         * (G[5]*G[6]) * Popu[i-1,5]

            # Coupling with rminusa (state 7)
            - k_rminusb_rminusa     * (G[6]*G[7]) * Popu[i-1,6]
            + k_rminusa_rminusb     * (G[7]*G[6]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 7: rminusa
        #-------------
        pump_rmin = (mode == :rmin) ? laser_here : 0
        Popu[i,7] = Popu[i-1,7] + 1/G[7] * (

            # Pumping if this mode is :rmin
            + pump_rmin

            # Relaxation
            - G[7] * relax_rminusa * Popu[i-1,7]

            # Coupling with rplus (state 1)
            - k_rminusa_rplus       * (G[7]*G[1]) * Popu[i-1,7]
            + k_rplus_rminusa       * (G[1]*G[7]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
            + k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusa_dfr         * (G[7]*G[3]) * Popu[i-1,7]
            + k_dfr_rminusa         * (G[3]*G[7]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusa_dminus      * (G[7]*G[4]) * Popu[i-1,7]
            + k_dminus_rminusa      * (G[4]*G[7]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusa_rfr        * (G[7]*G[5]) * Popu[i-1,7]
            + k_rfr_rminusa        * (G[5]*G[7]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rminusa_rminusb    * (G[7]*G[6]) * Popu[i-1,7]
            + k_rminusb_rminusa    * (G[6]*G[7]) * Popu[i-1,6]
        ) * dt

        #-------------
        # 7) Control
        # Checks if the balance of the populations is correct at each timestep
        #-------------

        # 1. Integral of Laser Input 
        control[i,1] = control[i-1,1] + laser_here * dt

        # 2. Overall Population of all substates
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )
        
        # 3. Relaxed Population of all substates
        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt

        # 4. Sum of Laser Input, Overall Population and Relaxed Population should equal to 0 at every timestep
        control[i,4] = control[i,1] - control[i,2] - control[i,3]
    end

    #---------------------------------------------------------
    # 8) Finally, compute the "Bleach" signal 
    #    for each state: (1 - 2*g[i]*Popu[i])^2
    #---------------------------------------------------------
    for k in 1:n
        Bleach[:,k] = (1 .- 2 .* G[k] .* Popu[:,k]) .^2
    end

# ------------------------------------------------------
    # 9) Dictionary (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1], t_0 = p[2])

    # B) Relaxation
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[3:9]
    )

    # C) Degeneracy
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G
    )

    # D) Couplings between states
    #    We differ between forward and backward couplings
    #    Syntax: (x, y, fwd_idx, bwd_idx)
    #      => t_x_y   => p[fwd_idx]
    #      => t_y_x   => p[bwd_idx]
    couplings = [
        (:rplus,       :dminusomega, 10, 10),
        (:rplus,       :dfr,         11, 11),
        (:rplus,       :dminus,      12, 12),
        (:rplus,       :rfr,         13, 28),  # asym
        (:rplus,       :rminusb,     14, 14),
        (:rplus,       :rminusa,     15, 15),

        (:dminusomega, :dfr,         16, 16),
        (:dminusomega, :dminus,      17, 17),
        (:dminusomega, :rfr,         18, 29),  # asym
        (:dminusomega, :rminusb,     19, 19),
        (:dminusomega, :rminusa,     20, 20),

        (:dfr,         :dminus,      21, 21),
        (:dfr,         :rfr,         22, 30),  # asym
        (:dfr,         :rminusb,     23, 23),
        (:dfr,         :rminusa,     24, 24),

        (:dminus,      :rfr,         25, 31),  # asym
        (:dminus,      :rminusb,     26, 26),
        (:dminus,      :rminusa,     27, 27),

        (:rfr,         :rminusb,     32, 34),  # asym
        (:rfr,         :rminusa,     33, 36),  # asym

        (:rminusb,     :rminusa,     35, 35)
    ]

    rows_all = Vector{Dict{Symbol,Any}}()
    for (x, y, fwd_idx, bwd_idx) in couplings
        push!(rows_all, Dict(
            :param_name => "t_$(x)_$(y)",
            :value      => p[fwd_idx],
            :ratio     => p[fwd_idx]/p[bwd_idx]

        ))
        push!(rows_all, Dict(
            :param_name => "t_$(y)_$(x)",
            :value      => p[bwd_idx],
            :ratio     => p[bwd_idx]/p[fwd_idx]
        ))
    end
    df_koppl_all = DataFrame(rows_all)

    # Function to create DataFrame for each mode
    function df_koppl_for_mode(m::Symbol)
        subrows = Vector{Dict{Symbol,Any}}()
        for (x,y, fwd_idx, bwd_idx) in couplings
            if x == m || y == m
                # x->y
                push!(subrows, Dict(
                    :param_name => "t_$(x)_$(y)",
                    :value      => p[fwd_idx],
                    :ratio     => p[fwd_idx]/p[bwd_idx]
                ))
                # y->x
                push!(subrows, Dict(
                    :param_name => "t_$(y)_$(x)",
                    :value      => p[bwd_idx],
                    :ratio     => p[bwd_idx]/p[fwd_idx]
                ))
            end
        end
        return DataFrame(subrows)
    end

    # Create DataFrame for each mode
    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Full output dictionary
    outputDict = Dict(
        :Model       => "Model 3.1",
        :fitted_mode => mode,
        :Parameter   => Dict(
            :All        => p,
            :Puls       => df_puls,
            :Relaxation => df_relax,
            :Entartung  => df_entartung,
            :Kopplung   => dict_koppl
        ),
        :Zeitschritt => dt
    )

    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end

"""
    The following function can fit the model 3.1 to the data. For more information on the model, see the model function model3_1.
    The function uses a global optimizer of the package "BlackBoxOptim.jl". It is recommend to use the :adaptive_de_rand_1_bin_radiuslimited algorithm 
    to have a robustfree fit. 
    This function fits three data sets to one set of parameters. The data sets are:
        - time_dmin, traces_dmin, use_dmin: data set for the dmin pumped experiment
        - time_rfr, traces_rfr, use_rfr: data set for the rfr pumped experiment
        - time_rmin, traces_rmin, use_rmin: data set for the rmin pumped experiment
    The time argument is the time vector of the corresponding data set, the traces argument is the data set itself 
    of the used traces and the use argument is a vector which contains the indices of the traces which should be used for the fit.
    The coha argument is an array. coha[1] contains the index of the mode where we observed the coherent artifact in the rfr pumped experiment.
    coha[2] contains the indices of the time_vector where the coherent artifact is observed.
    The function returns three ReservoirModel objects for the three data sets.
        julia> blackbox_model_3_1(
                time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
                time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
                time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
                coha;                           # cohrent artifact array
                dt=0.1                          # time step
               )
"""
function blackbox_model_3_1(
    time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
    time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
    time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
    coha;                           # cohrent artifact array
    dt=0.1                          # time step
)
    # 1) Define the search space
    lb = vcat(
        [0.0, -20.0], # α_dmin, t_0_dmin
        [0.0, -20.0], # α_rfr, t_0_rfr
        [0.0, -20.0], # α_rmin, t_0_rmin
        fill(1.0,7), # relaxation_params   
        fill(1.0,27), # kopplungs_params
        fill(1.0,7) # g_params
    )
    ub = vcat(
        [5.0, 20.0], # α_dmin, t_0_dmin
        [5.0, 20.0], # α_rfr, t_0_rfr
        [5.0, 20.0], # α_rmin, t_0_rmin
        fill(200000.0,7), # relaxation_params   
        fill(200000.0,27), # kopplungs_params
        2.0,              # g_rplus
        17.0,             # g_dminomega
        17.0,             # g_dfr
        17.0,             # g_dminus
        3.0,              # g_rfr
        5.0,              # g_rminusb
        2.0               # g_rminusa
    )
    # 2) Construct time vectors for each experiment
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    # 3) Define the cost function
    # The function "cost" computes a total cost based on the weighted squared differences
    # between the model predictions and experimental traces. It includes penalties to
    # constrain the simulated populations.
    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1];  t0_dmin  = p[2]
        α_rfr   = p[3];  t0_rfr   = p[4]
        α_rmin  = p[5];  t0_rmin  = p[6]

        relax_p = p[7:13]      # 7 relaxation
        koppl_p = p[14:40]     # 27 couplings
        g     = p[41:end]      # 7 degeneracy

        # 2) Construct the models
        # With the use argument we select a subset of the model to fit to the data

        # dmin pumped experiment
        model_dmin = model3_1(
            t_model_dmin,
            vcat(α_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # dmin subset
        model_dmin_used = model_dmin.Bleach[:, use_dmin]

        # rfr pumped experiment
        model_rfr = model3_1(
            t_model_rfr,
            vcat(α_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # rfr subset
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        
        # rmin pumped experiment
        model_rmin = model3_1(
            t_model_rmin,
            vcat(α_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # rmin subset
        model_rmin_used = model_rmin.Bleach[:, use_rmin]

        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
        #   Also, replace NaN and Inf values with a large number (1e12)
        pred_dmin = hcat(v_pred_dmin...)
        for i in 1:length(pred_dmin)
            if isinf(pred_dmin[i]) || isnan(pred_dmin[i])
                pred_dmin[i] = 1e12
            end
        end
        pred_rfr  = hcat(v_pred_rfr...)
        for i in 1:length(pred_rfr)
            if isinf(pred_rfr[i]) || isnan(pred_rfr[i])
                pred_rfr[i] = 1e12
            end
        end
        pred_rmin = hcat(v_pred_rmin...)
        for i in 1:length(pred_rmin)
            if isinf(pred_rmin[i]) || isnan(pred_rmin[i])
                pred_rmin[i] = 1e12
            end
        end
        # 5) Select time steps without coherence artifact in rfr pumped experiment
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights 
        # Use the time_weight function to create special weights for short delays, since we observe 
        # fast dynamics in the beginning of the experiment.
        w_dmin = ones(size(pred_dmin))
        w_pertrace_dmin = [
            5.0,  # rplus
            5.0,  # rfr
            1.0,  # rminusb
            6.0   # rminusa

        ]
        for m in 1:size(w_dmin,1)
            for n in 1:size(w_dmin,2)
                w_dmin[m,n] = time_weight(time_dmin[m],weight=3) * w_pertrace_dmin[n]
            end
        end
    
        w_rfr = ones(size(pred_rfr))
        w_pertrace_rfr = [
            5.0,  # rplus
            1.0,  # dfr
            1.0,  # rfr
            5.0   # rminusa
        ]
        for m in 1:size(w_rfr,1)
            for n in 1:size(w_rfr,2)
                if (n == coha[1]) && (m ∉ selected_idx_rfr)
                    w_rfr[m,n] = 0.0
                else
                    w_rfr[m,n] = time_weight(time_rfr[m],weight=3) * w_pertrace_rfr[n]
                end
            end
        end
        w_rmin = ones(size(pred_rmin))
        w_pertrace_rmin = [
            5.0,  # rplus
            1.0,  # dminusomega
            7.0,  # rfr
            1.0,  # rminusb
            5.0   # rminusa
        ]
        for m in 1:size(w_rmin,1)
            for n in 1:size(w_rmin,2)
                w_rmin[m,n] = time_weight(time_rmin[m],weight=3) * w_pertrace_rmin[n]
            end
        end
        
        # 7) Since we sum over all Data points, we need to normalize the weights
        # to the length of the given dataset
            length_dmin = length(time_dmin)
            length_rfr = length(time_rfr)
            length_rmin = length(time_rmin)
            norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
            norm_weight ./= maximum(norm_weight)
            cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
            cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
            cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 8) penalties
        # We define a maximum Population with max_P
        # If the Population is above this value, we add a penalty to the cost function
            max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
            P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
                norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
                norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
            total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty


        return total_cost
    end

    # Use bboptimize to find the best parameters
    println("Starting optimization...")

    result = bboptimize(
        cost; # cost function
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)], # search space
        Method       = :adaptive_de_rand_1_bin_radiuslimited ,   # or :cmaes, :genetic, etc.
        MaxSteps     = 10_000_000, # or some large budget
        PopulationSize = 2_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )
    println("Best solution = ", best_candidate(result))
    # 9) Extract the parameters
        ## Pumppulse
            ### dmin
                pulse_params_dmin  =  best_candidate(result)[1:2]
            ### rfr
                pulse_params_rfr  =  best_candidate(result)[3:4]
            ### rmin
                pulse_params_rmin  =  best_candidate(result)[5:6]
        ## Relaxation
            relaxation_params =  best_candidate(result)[7:13]
        ## Kopplung
            kopplungs_params = best_candidate(result)[14:40]
        ## Degeneracy
            g = best_candidate(result)[41:end]
    # 10) Create the model

    best_model_dmin  = model3_1(
        t_model_dmin,
        vcat(pulse_params_dmin,relaxation_params,kopplungs_params,g),
        mode=:dmin,
        dt=dt
    )
    best_model_rfr  = model3_1(
        t_model_rfr,
        vcat(pulse_params_rfr,relaxation_params,kopplungs_params,g),
        mode=:rfr,
        dt=dt
    )
    best_model_rmin  = model3_1(
        t_model_rmin,
        vcat(pulse_params_rmin,relaxation_params,kopplungs_params,g),
        mode=:rmin,
        dt=dt
    )

    return best_model_dmin,best_model_rfr,best_model_rmin
end

"""
    Modell 2.3
    Based on wavelength-scan measurements, it was observed that when pumping at "rmin", the system does not 
    selectively excite just one mode (rmina), but also simultaneously excites rminb. To account for this, an 
    additional pump pulse parameter was introduced. This parameter is always set to zero when "rmin" is not pumped.
    Model 2.3 is considered a refinement of Model 2.2.
    Modell 2.2:
    Model 2.2 describes the dynamics of all modes that exhibit continuous changes in experiments (ODT, UDT, HT).
    The modes included in this model are: rplus, dminusomega, dfr, dminus, rfr, rminusb, and rminusa.
    Each mode can be excited by a laser using the parameters α and t₀, which can be specified via the mode 
    keyword argument (e.g., mode=:dmin).
    Each mode can also undergo relaxation, controlled by parameters p[3:9], or interact with other modes 
    via coupling terms defined by p[10:30].
    Additionally, mode degeneracy can be set using parameters p[31:37].
"""
function model2_3(
    t,
    p::Vector{T};
    mode = :dmin,
    dt = 0.1
) where T<:Number

    #---------------------------------------------------------
    # 1) We have 7 excited states. We'll store their populations
    #    in Popu[:,1] through Popu[:,7].
    #    We'll also store a "bleach" signal in Bleach[:,1..7].
    #---------------------------------------------------------
    n = 7  # number of states
    Popu   = zeros(T, length(t), n)
    Bleach = zeros(T, length(t), n)
    control = zeros(T, length(t), 4)

    #---------------------------------------------------------
    # 2) Unpack parameters from p:
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    α_rb     = 0   # amplitude scaling of laser pulse for rminb
    if mode == :rmin
        α_rb = p[2]
    end
    t_0      = p[3]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width to observed pulse width of pump pulse

    # Relaxation rates to ground state 
    relax_rplus       = 1/p[4]
    relax_dminusomega = 1/p[5]
    relax_dfr         = 1/p[6]
    relax_dminus      = 1/p[7]
    relax_rfr         = 1/p[8]
    relax_rminusb     = 1/p[9]
    relax_rminusa     = 1/p[10]

    # Coupling constants between states 
    # For simplicity, I'll rename them 'k_x_y' to emphasize
    # they are coupling rates. Each will get multiplied by
    # the product of degeneracies g_x*g_y for symmetrical flow.

    k_rplus_dminusomega   = 1/p[11]
    k_rplus_dfr           = 1/p[12]
    k_rplus_dminus        = 1/p[13]
    k_rplus_rfr           = 1/p[14]
    k_rplus_rminusb       = 1/p[15]
    k_rplus_rminusa       = 1/p[16]

    k_dminusomega_rplus   = k_rplus_dminusomega
    k_dminusomega_dfr     = 1/p[17]
    k_dminusomega_dminus  = 1/p[18]
    k_dminusomega_rfr     = 1/p[19]
    k_dminusomega_rminusb = 1/p[20]
    k_dminusomega_rminusa = 1/p[21]

    k_dfr_rplus           = k_rplus_dfr
    k_dfr_dminusomega     = k_dminusomega_dfr
    k_dfr_dminus          = 1/p[22]
    k_dfr_rfr             = 1/p[23]
    k_dfr_rminusb         = 1/p[24]
    k_dfr_rminusa         = 1/p[25]

    k_dminus_rplus        = k_rplus_dminus
    k_dminus_dminusomega  = k_dminusomega_dminus
    k_dminus_dfr          = k_dfr_dminus
    k_dminus_rfr          = 1/p[26]
    k_dminus_rminusb      = 1/p[27]
    k_dminus_rminusa      = 1/p[28]

    k_rfr_rplus           = k_rplus_rfr       #1/p[29]
    k_rfr_dminusomega     = k_dminusomega_rfr #1/p[30]
    k_rfr_dfr             = k_dfr_rfr         #1/p[31]
    k_rfr_dminus          = k_dminus_rfr      #1/p[32]
    k_rfr_rminusb         = 1/p[29]
    k_rfr_rminusa         = 1/p[30]

    k_rminusb_rplus       = k_rplus_rminusb
    k_rminusb_dminusomega = k_dminusomega_rminusb
    k_rminusb_dfr         = k_dfr_rminusb
    k_rminusb_dminus      = k_dminus_rminusb
    k_rminusb_rfr         = k_rfr_rminusb    #1/p[35]
    k_rminusb_rminusa     = 1/p[31]

    k_rminusa_rplus       = k_rplus_rminusa
    k_rminusa_dminusomega = k_dminusomega_rminusa
    k_rminusa_dfr         = k_dfr_rminusa
    k_rminusa_dminus      = k_dminus_rminusa
    k_rminusa_rfr         = k_rfr_rminusa    #1/p[37]
    k_rminusa_rminusb     = k_rminusb_rminusa

    # 3) Degeneracy of the states:
    # only fixed values for rplus and rfr
    g_rplus   = 1.0 #p[32]
    g_d_omega = p[33]
    g_dfr     = p[34]
    g_dminus  = p[35]
    g_rfr     = 1.0 #p[36]
    g_rminb   = p[37]
    g_rmina   = p[38]

    # We'll store them in an array so we can do bleach:
    # State order: 
    #   1->rplus, 2->dminusomega, 3->dfr, 4->dminus, 
    #   5->rfr,   6->rminusb,     7->rminusa
    G = [g_rplus, g_d_omega, g_dfr, g_dminus, g_rfr, g_rminb, g_rmina]

    #---------------------------------------------------------
    # 4) Precompute the pulse shape at each time point
    f_g    = gauss.(t,Ref([σ,t_0]))

    #---------------------------------------------------------
    # 5) Mode-specific initial conditions:
    #    :dmin, :rmin, or :rfr  -> which state is initially pumped?
    #---------------------------------------------------------
        pump_first = α * f_g[1] * dt
    # Initially, one state is pumped:
        Popu[1,1] = (mode == :rplus)     ? 1/G[1] * pump_first         : 0
        Popu[1,2] = (mode == :dminomega) ? 1/G[2] * pump_first         : 0
        Popu[1,3] = (mode == :dfr)       ? 1/G[3] * pump_first         : 0
        Popu[1,4] = (mode == :dmin)      ? 1/G[4] * pump_first         : 0
        Popu[1,5] = (mode == :rfr)       ? 1/G[5] * pump_first         : 0
        Popu[1,6] = (mode == :rmin)      ? 1/G[6] * α_rb * f_g[1] * dt : 0
        Popu[1,7] = (mode == :rmin)      ? 1/G[7] * pump_first         : 0
        control[1,1] = pump_first + α_rb * f_g[1] * dt 
        control[1,2] = G[1] * Popu[1,1] + (
            + G[2] * Popu[1,2] 
            + G[3] * Popu[1,3] 
            + G[4] * Popu[1,4] 
            + G[5] * Popu[1,5] 
            + G[6] * Popu[1,6] 
            + G[7] * Popu[1,7]
        )


    #---------------------------------------------------------
    # 6) Time stepping with explicit Euler:
    #    We'll update all 7 states at each time step.
    #
    #    Key change for symmetrical coupling:
    #      * For a coupling k_x_y between state x and y,
    #        we multiply by (G[x]*G[y]) in BOTH equations:
    #
    #        Popu[i,x] += + k_xy * G[x]*G[y]*(Popu[i-1,y] - Popu[i-1,x])
    #        Popu[i,y] += + k_xy * G[x]*G[y]*(Popu[i-1,x] - Popu[i-1,y])
    #
    #    The same factor appears in both directions.
    #---------------------------------------------------------
    for i in 2:length(t)

        # Precompute the laser amplitude at this time, if needed
        laser_here = α * f_g[i]

        #-------------
        # State 1: rplus 
        #-------------
        # If this is the pumped state in :rplus mode, we add laser_here
        # (same code runs for all modes, but laser_here=0 if not used)
        pump_rplus = (mode == :rplus) ?  laser_here : 0
        Popu[i,1] = Popu[i-1,1] + 1/G[1] * (

            # Pumping if this mode is :rplus
            + pump_rplus

            # Relaxation to ground
            - G[1] * relax_rplus * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rplus_dminusomega * (G[1]*G[2]) * Popu[i-1,1]
            + k_dminusomega_rplus * (G[2]*G[1]) * Popu[i-1,2] 

            # Coupling with dfr (state 3)
            - k_rplus_dfr         * (G[1]*G[3]) * Popu[i-1,1]
            + k_dfr_rplus         * (G[3]*G[1]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rplus_dminus      * (G[1]*G[4]) * Popu[i-1,1]
            + k_dminus_rplus      * (G[4]*G[1]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rplus_rfr         * (G[1]*G[5]) * Popu[i-1,1]
            + k_rfr_rplus         * (G[5]*G[1]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rplus_rminusb     * (G[1]*G[6]) * Popu[i-1,1]
            + k_rminusb_rplus     * (G[6]*G[1]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rplus_rminusa     * (G[1]*G[7]) * Popu[i-1,1]
            + k_rminusa_rplus     * (G[7]*G[1]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 2: dminusomega
        #-------------
        pump_dminomega = (mode == :dminomega) ? laser_here : 0
        Popu[i,2] = Popu[i-1,2] + 1/G[2] * (

            # Pumping if this mode is :dminomega
            + pump_dminomega

            # Relaxation
            - G[2] * relax_dminusomega * Popu[i-1,2]

            # Coupling with rplus (state 1)
            - k_dminusomega_rplus   * (G[2]*G[1]) * Popu[i-1,2]
            + k_rplus_dminusomega   * (G[1]*G[2]) * Popu[i-1,1]

            # Coupling with dfr (state 3)
            - k_dminusomega_dfr     * (G[2]*G[3]) * Popu[i-1,2]
            + k_dfr_dminusomega     * (G[3]*G[2]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_dminusomega_dminus  * (G[2]*G[4]) * Popu[i-1,2]
            + k_dminus_dminusomega  * (G[4]*G[2]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dminusomega_rfr     * (G[2]*G[5]) * Popu[i-1,2]
            + k_rfr_dminusomega     * (G[5]*G[2]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]
            + k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]
            + k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 3: dfr
        #-------------
        pump_dfr = (mode == :dfr) ? laser_here : 0
        Popu[i,3] = Popu[i-1,3] + 1/G[3] * (

            # Pumping if this mode is :dfr
            + pump_dfr

            # Relaxation
            - G[3] * relax_dfr * Popu[i-1,3]

            # Coupling with rplus (state 1)
            - k_dfr_rplus       * (G[3]*G[1]) * Popu[i-1,3]
            + k_rplus_dfr       * (G[1]*G[3]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dfr_dminusomega * (G[3]*G[2]) * Popu[i-1,3]
            + k_dminusomega_dfr * (G[2]*G[3]) * Popu[i-1,2]

            # Coupling with dminus (state 4)
            - k_dfr_dminus      * (G[3]*G[4]) * Popu[i-1,3]
            + k_dminus_dfr      * (G[4]*G[3]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]
            + k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dfr_rminusb     * (G[3]*G[6]) * Popu[i-1,3]
            + k_rminusb_dfr     * (G[6]*G[3]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dfr_rminusa     * (G[3]*G[7]) * Popu[i-1,3]
            + k_rminusa_dfr     * (G[7]*G[3]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 4: dminus
        #-------------
        pump_dmin = (mode == :dmin) ? laser_here : 0
        Popu[i,4] = Popu[i-1,4] + 1/G[4] * (

            # Pumping if this mode is :dmin
            + pump_dmin

            # Relaxation
            - G[4] * relax_dminus * Popu[i-1,4]

            # Coupling with rplus (state 1)
            - k_dminus_rplus       * (G[4]*G[1]) * Popu[i-1,4]
            + k_rplus_dminus       * (G[1]*G[4]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_dminus_dminusomega * (G[4]*G[2]) * Popu[i-1,4]
            + k_dminusomega_dminus * (G[2]*G[4]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_dminus_dfr         * (G[4]*G[3]) * Popu[i-1,4]
            + k_dfr_dminus         * (G[3]*G[4]) * Popu[i-1,3]

            # Coupling with rfr (state 5)
            - k_dminus_rfr         * (G[4]*G[5]) * Popu[i-1,4]
            + k_rfr_dminus         * (G[5]*G[4]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_dminus_rminusb     * (G[4]*G[6]) * Popu[i-1,4]
            + k_rminusb_dminus     * (G[6]*G[4]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_dminus_rminusa     * (G[4]*G[7]) * Popu[i-1,4]
            + k_rminusa_dminus     * (G[7]*G[4]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 5: rfr
        #-------------
        pump_rfr = (mode == :rfr) ? laser_here : 0
        Popu[i,5] = Popu[i-1,5] + 1/G[5] * (
            # Pumping if this mode is :rfr
            + pump_rfr
            
            # Relaxation
            - G[5] * relax_rfr * Popu[i-1,5]

            # Coupling with rplus (state 1)
            - k_rfr_rplus       * (G[5]*G[1]) * Popu[i-1,5]
            + k_rplus_rfr       * (G[1]*G[5]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rfr_dminusomega * (G[5]*G[2]) * Popu[i-1,5]
            + k_dminusomega_rfr * (G[2]*G[5]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rfr_dfr         * (G[5]*G[3]) * Popu[i-1,5]
            + k_dfr_rfr         * (G[3]*G[5]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rfr_dminus      * (G[5]*G[4]) * Popu[i-1,5]
            + k_dminus_rfr      * (G[4]*G[5]) * Popu[i-1,4]

            # Coupling with rminusb (state 6)
            - k_rfr_rminusb     * (G[5]*G[6]) * Popu[i-1,5]
            + k_rminusb_rfr     * (G[6]*G[5]) * Popu[i-1,6]

            # Coupling with rminusa (state 7)
            - k_rfr_rminusa     * (G[5]*G[7]) * Popu[i-1,5]
            + k_rminusa_rfr     * (G[7]*G[5]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 6: rminusb
        #-------------
        #pump_rminb = (mode == :rminb) ? laser_here : 0
        pump_rminb  = α_rb * f_g[i] * dt
        Popu[i,6] = Popu[i-1,6] + 1/G[6] * (

            # Pumping if this mode is :rminb
            + pump_rminb

            # Relaxation
            - G[6] * relax_rminusb * Popu[i-1,6]

            # Coupling with rplus (state 1)
            - k_rminusb_rplus       * (G[6]*G[1]) * Popu[i-1,6]
            + k_rplus_rminusb       * (G[1]*G[6]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusb_dminusomega * (G[6]*G[2]) * Popu[i-1,6]
            + k_dminusomega_rminusb * (G[2]*G[6]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusb_dfr         * (G[6]*G[3]) * Popu[i-1,6]
            + k_dfr_rminusb         * (G[3]*G[6]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusb_dminus      * (G[6]*G[4]) * Popu[i-1,6]
            + k_dminus_rminusb      * (G[4]*G[6]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusb_rfr         * (G[6]*G[5]) * Popu[i-1,6]
            + k_rfr_rminusb         * (G[5]*G[6]) * Popu[i-1,5]

            # Coupling with rminusa (state 7)
            - k_rminusb_rminusa     * (G[6]*G[7]) * Popu[i-1,6]
            + k_rminusa_rminusb     * (G[7]*G[6]) * Popu[i-1,7]
        ) * dt

        #-------------
        # State 7: rminusa
        #-------------
        pump_rmin = (mode == :rmin) ? laser_here : 0
        Popu[i,7] = Popu[i-1,7] + 1/G[7] * (

            # Pumping if this mode is :rmin
            + pump_rmin

            # Relaxation
            - G[7] * relax_rminusa * Popu[i-1,7]

            # Coupling with rplus (state 1)
            - k_rminusa_rplus       * (G[7]*G[1]) * Popu[i-1,7]
            + k_rplus_rminusa       * (G[1]*G[7]) * Popu[i-1,1]

            # Coupling with dminusomega (state 2)
            - k_rminusa_dminusomega * (G[7]*G[2]) * Popu[i-1,7]
            + k_dminusomega_rminusa * (G[2]*G[7]) * Popu[i-1,2]

            # Coupling with dfr (state 3)
            - k_rminusa_dfr         * (G[7]*G[3]) * Popu[i-1,7]
            + k_dfr_rminusa         * (G[3]*G[7]) * Popu[i-1,3]

            # Coupling with dminus (state 4)
            - k_rminusa_dminus      * (G[7]*G[4]) * Popu[i-1,7]
            + k_dminus_rminusa      * (G[4]*G[7]) * Popu[i-1,4]

            # Coupling with rfr (state 5)
            - k_rminusa_rfr        * (G[7]*G[5]) * Popu[i-1,7]
            + k_rfr_rminusa        * (G[5]*G[7]) * Popu[i-1,5]

            # Coupling with rminusb (state 6)
            - k_rminusa_rminusb    * (G[7]*G[6]) * Popu[i-1,7]
            + k_rminusb_rminusa    * (G[6]*G[7]) * Popu[i-1,6]
        ) * dt

        #-------------
        # 7) Control
        # Checks if the balance of the populations is correct at each timestep
        #-------------
        
        # 1. Integral of Laser Input 
        control[i,1] = control[i-1,1] + (laser_here + α_rb * f_g[i]) * dt

        # 2. Overall Population of all substates 
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )

        # 3. Relaxed Population of all substates
        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt

        # 4. Sum of Laser Input, Overall Population and Relaxed Population should equal to 0 at every timestep
        control[i,4] = control[i,1] - control[i,2] - control[i,3]
    end

    #---------------------------------------------------------
    # 8) Finally, compute the "Bleach" signal 
    #    for each state: (1 - 2*g[i]*Popu[i])^2
    #---------------------------------------------------------
    for k in 1:n
        Bleach[:,k] = (1 .- 2 .* G[k] .* Popu[:,k]) .^2
    end

# ------------------------------------------------------
    # 9) Dictionary (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1],alpha_rb = p[2], t_0 = p[3])

    # B) Relaxation
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[4:10]
    )

    # C) Entartung
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G
    )

    # D) Couplings between states
    #    We differ between forward and backward couplings even though they are the same.
    #    Syntax: (x, y, fwd_idx, bwd_idx)
    #      => t_x_y   => p[fwd_idx]
    #      => t_y_x   => p[bwd_idx]
    couplings = [
        (:rplus,       :dminusomega, 11, 11),
        (:rplus,       :dfr,         12, 12),
        (:rplus,       :dminus,      13, 13),
        (:rplus,       :rfr,         14, 14),  
        (:rplus,       :rminusb,     15, 15),
        (:rplus,       :rminusa,     16, 16),

        (:dminusomega, :dfr,         17, 17),
        (:dminusomega, :dminus,      18, 18),
        (:dminusomega, :rfr,         19, 19),  
        (:dminusomega, :rminusb,     20, 20),
        (:dminusomega, :rminusa,     21, 21),

        (:dfr,         :dminus,      22, 22),
        (:dfr,         :rfr,         23, 23),  
        (:dfr,         :rminusb,     24, 24),
        (:dfr,         :rminusa,     25, 25),

        (:dminus,      :rfr,         26, 26),  
        (:dminus,      :rminusb,     27, 27),
        (:dminus,      :rminusa,     28, 28),

        (:rfr,         :rminusb,     29, 29),  
        (:rfr,         :rminusa,     30, 30),  

        (:rminusb,     :rminusa,     31, 31)
    ]

    rows_all = Vector{Dict{Symbol,Any}}()
    for (x, y, fwd_idx, bwd_idx) in couplings
        push!(rows_all, Dict(
            :param_name => "t_$(x)_$(y)",
            :value      => p[fwd_idx],
            :ratio     => p[fwd_idx]/p[bwd_idx]

        ))
        push!(rows_all, Dict(
            :param_name => "t_$(y)_$(x)",
            :value      => p[bwd_idx],
            :ratio     => p[bwd_idx]/p[fwd_idx]
        ))
    end
    df_koppl_all = DataFrame(rows_all)

    # Function to create DataFrame for each mode
    function df_koppl_for_mode(m::Symbol)
        subrows = Vector{Dict{Symbol,Any}}()
        for (x,y, fwd_idx, bwd_idx) in couplings
            if x == m || y == m
                # x->y
                push!(subrows, Dict(
                    :param_name => "t_$(x)_$(y)",
                    :value      => p[fwd_idx],
                    :ratio     => p[fwd_idx]/p[bwd_idx]
                ))
                # y->x
                push!(subrows, Dict(
                    :param_name => "t_$(y)_$(x)",
                    :value      => p[bwd_idx],
                    :ratio     => p[bwd_idx]/p[fwd_idx]
                ))
            end
        end
        return DataFrame(subrows)
    end
    # Create DataFrame for each mode
    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Full output dictionary
    outputDict = Dict(
        :Model       => "Model 2.3",
        :fitted_mode => mode,
        :Parameter   => Dict(
            :All        => p,
            :Puls       => df_puls,
            :Relaxation => df_relax,
            :Entartung  => df_entartung,
            :Kopplung   => dict_koppl
        ),
        :Zeitschritt => dt
    )

    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end

"""
    The following function can fit the model 2.3 to the data. For more information on the model, see the model function model2_3.
    The function uses a global optimizer of the package "BlackBoxOptim.jl". It is recommend to use the :adaptive_de_rand_1_bin_radiuslimited algorithm 
    to have a robustfree fit. 
    This function fits three data sets to one set of parameters. The data sets are:
        - time_dmin, traces_dmin, use_dmin: data set for the dmin pumped experiment
        - time_rfr, traces_rfr, use_rfr: data set for the rfr pumped experiment
        - time_rmin, traces_rmin, use_rmin: data set for the rmin pumped experiment
    The time argument is the time vector of the corresponding data set, the traces argument is the data set itself 
    of the used traces and the use argument is a vector which contains the indices of the traces which should be used for the fit.
    The coha argument is an array. coha[1] contains the index of the mode where we observed the coherent artifact in the rfr pumped experiment.
    coha[2] contains the indices of the time_vector where the coherent artifact is observed.
    The function returns three ReservoirModel objects for the three data sets.
        julia> blackbox_model_3_2(
                time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
                time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
                time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
                coha;                           # cohrent artifact array
                dt=0.1                          # time step
               )
"""
function blackbox_model_2_3(
    time_dmin,traces_dmin,use_dmin, # data set for the dmin pumped experiment
    time_rfr,traces_rfr,use_rfr,    # data set for the rfr pumped experiment
    time_rmin,traces_rmin,use_rmin, # data set for the rmin pumped experiment
    coha;                           # cohrent artifact array
    dt=0.1                          # time step
)
    # 1) Define the search space
    lb = vcat(
        [0.0, 0.0, -20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [0.0, 0.0, -20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [0.0, 0.0, -20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(1.0,7), # relaxation_params   
        fill(1.0,21), # kopplungs_params
        1.0,              # g_rplus
        1.0,             # g_dminomega
        1.0,             # g_dfr
        16.0,             # g_dminus
        1.0,              # g_rfr
        1.0,              # g_rminusb
        1.0               # g_rminusa
    )
    ub = vcat(
        [5.0, 5.0, 20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [5.0, 5.0, 20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [5.0, 5.0, 20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(20000.0,7), # relaxation_params   
        fill(20000.0,21), # kopplungs_params
        5.0,              # g_rplus
        17.0,             # g_dminomega
        17.0,             # g_dfr
        17.0,             # g_dminus
        5.0,              # g_rfr
        17.0,              # g_rminusb
        2.0               # g_rminusa
    )
    # 2) Construct time vectors for each experiment
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    # 3) Define the cost function
    # The function "cost" computes a total cost based on the weighted squared differences
    # between the model predictions and experimental traces. It includes penalties to
    # constrain the simulated populations.
    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1]; α_rb_dmin  = p[2];  t0_dmin  = p[3]
        α_rfr   = p[4]; α_rb_rfr   = p[5];  t0_rfr   = p[6]
        α_rmin  = p[7]; α_rb_rmin  = p[8];  t0_rmin  = p[9]

        relax_p = p[10:16]      # 7 relaxation
        koppl_p = p[17:37]     # 21 couplings
        g     = p[38:end]      # 7 degeneracy

        # 2) Construct the model
        # With the use argument we select a subset of the model to fit to the data

        # dmin pumped experiment
        model_dmin = model2_3(
            t_model_dmin,
            vcat(α_dmin,α_rb_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # dmin subset
        model_dmin_used = model_dmin.Bleach[:, use_dmin]

        # rfr pumped experiment
        model_rfr = model2_3(
            t_model_rfr,
            vcat(α_rfr, α_rb_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # rfr subset
        model_rfr_used = model_rfr.Bleach[:, use_rfr]

        # rmin pumped experiment
        model_rmin = model2_3(
            t_model_rmin,
            vcat(α_rmin,α_rb_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # rmin subset
        model_rmin_used = model_rmin.Bleach[:, use_rmin]

        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
        #   Also, replace NaN and Inf values with a large number (1e12)
        pred_dmin = hcat(v_pred_dmin...)
        for i in 1:length(pred_dmin)
            if isinf(pred_dmin[i]) || isnan(pred_dmin[i])
                pred_dmin[i] = 1e12
            end
        end
        pred_rfr  = hcat(v_pred_rfr...)
        for i in 1:length(pred_rfr)
            if isinf(pred_rfr[i]) || isnan(pred_rfr[i])
                pred_rfr[i] = 1e12
            end
        end
        pred_rmin = hcat(v_pred_rmin...)
        for i in 1:length(pred_rmin)
            if isinf(pred_rmin[i]) || isnan(pred_rmin[i])
                pred_rmin[i] = 1e12
            end
        end
        # 5) Select time steps without coherence artifact in rfr pumped experiment
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights 
        # Use the time_weight function to create special weights for short delays, since we observe 
        # fast dynamics in the beginning of the experiment.
        w_dmin = ones(size(pred_dmin))
        w_pertrace_dmin = [
            5.0,  # rplus
            5.0,  # rfr
            2.0,  # rminusb
            6.0   # rminusa

        ]
        for m in 1:size(w_dmin,1)
            for n in 1:size(w_dmin,2)
                w_dmin[m,n] = time_weight(time_dmin[m],weight=3) * w_pertrace_dmin[n]
            end
        end
    
        w_rfr = ones(size(pred_rfr))
        w_pertrace_rfr = [
            5.0,  # rplus
            1.0,  # dfr
            1.0,  # rfr
            5.0   # rminusa
        ]
        for m in 1:size(w_rfr,1)
            for n in 1:size(w_rfr,2)
                if (n == coha[1]) && (m ∉ selected_idx_rfr)
                    w_rfr[m,n] = 0.0
                else
                    w_rfr[m,n] = time_weight(time_rfr[m],weight=3) * w_pertrace_rfr[n]
                end
            end
        end
        w_rmin = ones(size(pred_rmin))
        w_pertrace_rmin = [
            5.0,  # rplus
            1.0,  # dminusomega
            7.0,  # rfr
            2.0,  # rminusb
            5.0   # rminusa
        ]
        for m in 1:size(w_rmin,1)
            for n in 1:size(w_rmin,2)
                w_rmin[m,n] = time_weight(time_rmin[m],weight=3) * w_pertrace_rmin[n]
            end
        end
        
        
        # 7) Since we sum over all Data points, we need to normalize the weights
        # to the length of the given dataset
            length_dmin = length(time_dmin)
            length_rfr = length(time_rfr)
            length_rmin = length(time_rmin)
            norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
            norm_weight ./= maximum(norm_weight)
            cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
            cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
            cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 8) Penaltys
        # We define a maximum Population with max_P
        # If the Population is above this value, we add a penalty to the cost function
            max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
            P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
                norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
                norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
            total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty

        return total_cost
    end

    # Use bboptimize to find the best parameters
    println("Starting optimization...")

    result = bboptimize(
        cost; # cost function
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)], #search space
        Method       = :dxnes,#:adaptive_de_rand_1_bin_radiuslimited ,   # or :cmaes, :genetic, etc.
        MaxSteps     = 100_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )
    println("Best solution = ", best_candidate(result))
    # 9) Extract the parameters
        ## Pumppulse
            ### dmin
                pulse_params_dmin  =  best_candidate(result)[1:3]
            ### rfr
                pulse_params_rfr  =  best_candidate(result)[4:6]
            ### rmin
                pulse_params_rmin  =  best_candidate(result)[7:9]
        ## Relaxation
            relaxation_params =  best_candidate(result)[10:16]
        ## Kopplung
            kopplungs_params = best_candidate(result)[17:37]
        ## Degeneracy
            g = best_candidate(result)[38:end]

    # 10) Create the the fitted model

    best_model_dmin  = model2_3(
        t_model_dmin,
        vcat(pulse_params_dmin,relaxation_params,kopplungs_params,g),
        mode=:dmin,
        dt=dt
    )
    best_model_rfr  = model2_3(
        t_model_rfr,
        vcat(pulse_params_rfr,relaxation_params,kopplungs_params,g),
        mode=:rfr,
        dt=dt
    )
    best_model_rmin  = model2_3(
        t_model_rmin,
        vcat(pulse_params_rmin,relaxation_params,kopplungs_params,g),
        mode=:rmin,
        dt=dt
    )

    return best_model_dmin,best_model_rfr,best_model_rmin
end

"""
This Funktion simulates the relaxation of a vibrational mode into their groundstate after laser excitation by a gaussian pulse.
dn/dt = α * gauss(t) - k * n(t) \\
    -> n(t)  = n(t-1) +  (α * gauss(t) - k * n(t-1)) * dt \\ 
function model1(t,p::Vector{T}; dt=0.1) where T<:Number

    α,k,t_0,offset = p 
    σ = 6.67
    n = zeros(T,length(t))
    f_g = gauss.(t,Ref([σ,t_0]))
    for i in eachindex(t)
        if i == 1
            n[1] = dt*(α * f_g[1])
        else 
            n[i] = n[i-1]+ (α * f_g[i] - 1/k * n[i-1]) * dt
        end
    end
    bleach = (1 .- n).^2
    return  bleach
end

"""
function model1(t,p::Vector{T}; dt=0.1) where T<:Number

    α,t_0,k = p
    σ = 6.67
    n = zeros(T,length(t))
    f_g = gauss.(t,Ref([σ,t_0]))
    for i in eachindex(t)
        if i == 1
            n[1] = dt*(α * f_g[1])
        else 
            n[i] = n[i-1]+ (α * f_g[i] - k * n[i-1]) * dt
        end
    end
    bleach = (1 .- n).^2
    return  ReservoirModel(n,bleach,p,t,1,1,Dict())
end


"""
    This function is used to create a time dependent weight for the cost function.
        julia> time_weigtht(3; t1=0, t2=5, weight=3.0)
        3.0
        julia> time_weigtht(6; t1=0, t2=5, weight=3.0)
        1.0
"""
function time_weight(t;
    t1 = -15,
    t2 = 30,
    weight = 3.0
    )
    if t1 <= t < t2
        return weight  # heavier weight for early times
    else
        return 1.0   # normal weight for later times
    end
end


"""
Gauss Function: \\
function gauss(p,t) \\   
    σ,t_0 = p \\
    1/(σ*sqrt(2π)) * exp(-(t-t_0)^2/(2σ^2)) \\
julia > t = collect(-10:0.1:10) \\
julia > f_g = [gauss([1,0],_t) for _t in t] 
"""
function gauss(t,p::AbstractVector{T})  where T
    σ,t_0 = p
    
	return 1.0 / (σ * sqrt(2π)) * exp.(-(t - t_0)^2 / (2σ^2))
end
