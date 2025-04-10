"""
    Model 3.3
    Modell 3.3 ist eine Weiterentwicklung von Modell 3.2. Es zeigt als wären gewissen Raten so langsam, dass sie aus der Simulation entfernt werden können.
    unter model.Dictionary sind die entfernte, bzw. festgelegten Parameter aufgelistet.
    Modell 3.2
    Aus dem Wellenlängenscan geht hervor, dass wenn wir "rmin" pumpen, wir nicht diskret rmina pumpen, sondern auch gleichzeitig rminb. Weswegen hier ein zusätzlicher 
    Pumppulsparameter hinzugefügt wird, elcher immer 0 ist wenn rmin nicht gepumpt wird.
    Wir betrachten Modell 3.2 als Weiterentwicklung von Modell 3.1.
    Modell 3.1:
    Das Modell 3.1 ist eine Weiterentwicklung des Modells 2.2. Modell 2.2 versagt am deutlichsten bei der Beschreibung
    der rfr Mode. Hier kann keine schnelle Dynamik erreicht werden. Das Modell 3.1 versucht dies zu verbessern, indem
    es die Hin und R"uckrate zwischen der rfr Mode nicht mehr als gleich annimmt.
    Wie in Modell 2.2 beschreiben die Parameter p[1:2] die Laseranregung. Die Parameter p[3:9] beschreiben die Relaxationsraten
    der Moden. Die Parameter p[10:50] beschreiben die Kopplungsraten zwischen den Moden. Die Parameter p[51:57] beschreiben die
    Entartung der Moden.
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
    #    p[1], p[2], etc.  (Keep same indexing as your code.)
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    α_rb     = 0   # amplitude scaling of laser pulse for rminb
    if mode == :rmin
        α_rb = p[2]
    end
    t_0      = p[3]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width (was hard-coded)

    # Relaxation rates (to ground state) 
    relax_rplus       = 1/p[4]
    #relax_dminusomega = 1/p[5]
    #relax_dfr         = 1/p[6]
    relax_dminus      = 1/p[5]
    #relax_rfr         = 1/p[4]
    relax_rminusb     = 1/p[6]
    relax_rminusa     = 1/p[7]

    # Coupling constants between states (original indexing)
    # For simplicity, I'll rename them 'k_xy' to emphasize
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
    # p[31] is g_d, while g_r, g_d_omega are hard-coded to 1 below
        g_rplus   = 1 #p[38]
        g_d_omega = 2
        g_dfr     = 12
        g_dminus  = 16
        g_rfr     = 1
        g_rminb   = 3
        g_rmina   = 2 #p[44]


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
    # 7) Time stepping with explicit Euler:
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
        # Control
        #-------------

        control[i,1] = control[i-1,1] + (laser_here + α_rb * f_g[i]) * dt
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )

        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            #+ G[2] * relax_dminusomega * Popu[i-1,2]
            #+ G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            #+ G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        
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
    # 9) Dictionary aufbauen (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1],alpha_rb = p[2], t_0 = p[3])

    # B) Relaxation
    #    p[3..9], 7 Parameter => rplus, dminusomega, dfr, dminus, rfr, rminusb, rminusa
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

    # C) Entartung
    #    p[37..43] => 7 Parameter
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G#p[38:44]
    )

    # D) Kopplungen (jetzt mit asymmetrischen Einträgen)
    #    Wir müssen nun Vorwärts- und Rückwärts-Parameter explizit angeben.
    #    Syntax: (x, y, fwd_idx, bwd_idx)
    #      => t_x_y   => p[fwd_idx]
    #      => t_y_x   => p[bwd_idx]
    couplings = [
        (:rplus,       :dminusomega, 8, 8),
        (:rplus,       :dfr,         9, 9),
        (:rplus,       :rminusb,     10, 10), # asyn

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

    # pro Mode eine gefilterte Kopplungstabelle
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

    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end
    # E) Entfernte Parameter
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
    # Gesamtdictionary
    outputDict = Dict(
        :Model       => "Model 3.3 - ODT",
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

    # ------------------------------------------------------
    # 10) Rückgabe als ReservoirModel (inkl. Dictionary)
    # ------------------------------------------------------
    # Falls dein ReservoirModel-Struct kein Dictionary-Feld hat,
    # kannst du hier natürlich nur das Modell und das Dict getrennt zurückgeben.
    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end


function blackbox_model_3_3(
    time_dmin,traces_dmin,use_dmin,
    time_rfr,traces_rfr,use_rfr,
    time_rmin,traces_rmin,use_rmin,
    coha;
    dt=0.1
)
    # 1) Define priors
    initial_params = [    
        1.1782,  # α_dmin
        7.196e-7,# α_rb_dmin
        -6.4929, # t_0_dmin
        0.207,   # α_rfr
        9.8e-7,  # α_rb_rfr
        -5.92,   # t_0_rfr
        0.0389,  # α_rmin
        2.8652  ,# α_rb_rmin
        -5.922,  # t_0_rmin
        3759.4,  # relax_rplus
        7.9595,  # relax_dminus
        65.375,  # relax_rminusb
        25.293,  # relax_rminusa
        16.302,  # rplus_dminusomega
        50.119,  # rplus_dfr
        28.914,  # rplus_rminusb
        5000.0, # dminusomega_rfr
        295.47,  # dminusomega_rminusa
        122.49,  # dfr_rminusb
        24.151,  # dminus_rfr
        3.1292,  # rfr_dminusomega
        5000.0, # rfr_dminus
        5000.0, # rfr_rminusb
        20.925,  # rfr_rminusa
        11.61,   # rminusb_rfr
        5000.0, # rminusa_rfr
    ]
    # 2) Define the search space
    lb = vcat(
        [0.0, 0.0, -20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [0.0, 0.0, -20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [0.0, 0.0, -20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(1.0,4), # relaxation_params   
        fill(1.0,13), # kopplungs_params
    )
    ub = vcat(
        [50.0, 1e-6, 20.0], # α_dmin, α_rb_dmin, t_0_dmin
        [50.0, 1e-6, 20.0], # α_rfr,  α_rb_rfr,  t_0_rfr
        [50.0, 50.0, 20.0], # α_rmin, α_rb_rmin, t_0_rmin
        fill(5000.0,4), # relaxation_params   
        fill(5000.0,13), # kopplungs_params
    )
    # 3) Construct time points for each mode
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1]; α_rb_dmin  = p[2];  t0_dmin  = p[3]
        α_rfr   = p[4]; α_rb_rfr   = p[5];  t0_rfr   = p[6]
        α_rmin  = p[7]; α_rb_rmin  = p[8];  t0_rmin  = p[9]

        relax_p = p[10:13]     # 4 relaxation
        koppl_p = p[14:26]     # 13 couplings
        g     = 1#p[28:end]    # 4 degeneracy

        # 2) Construct the model
        model_dmin = model3_3(
            t_model_dmin,
            vcat(α_dmin,α_rb_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_dmin_used = model_dmin.Bleach[:, use_dmin]
        model_rfr = model3_3(
            t_model_rfr,
            vcat(α_rfr, α_rb_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        model_rmin = model3_3(
            t_model_rmin,
            vcat(α_rmin,α_rb_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rmin_used = model_rmin.Bleach[:, use_rmin]
        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
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
        # 5) Select time steps without coherence artifact in rfr
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights
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
        
        # 6) Data likelihood
            length_dmin = length(time_dmin)
            length_rfr = length(time_rfr)
            length_rmin = length(time_rmin)
            norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
            norm_weight ./= maximum(norm_weight)
            cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
            cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
            cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 7) Penaltys
            max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
            P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
                norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
                norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
            total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty

        return total_cost
    end
    println("Starting optimization...")

    result = bboptimize(
        cost, # cost functionw_
        initial_params;          # initial guess
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)],
        Method       = :dxnes,#:adaptive_de_rand_1_bin_radiuslimited , #dxnes, or :cmaes, :genetic, etc.
        #NThreads=Threads.nthreads()-1,
        MaxSteps     = 1_000_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )

    #println("Finished. Best fitness = ", minimum_fitness(result))
    println("Best solution = ", best_candidate(result))
    # 7) Extract the parameters
        ## Pumppulse
            ### dmin
                pulse_params_dmin  =  best_candidate(result)[1:3]
            ### rfr
                pulse_params_rfr  =  best_candidate(result)[4:6]
            ### rmin
                pulse_params_rmin  =  best_candidate(result)[7:9]
        ## Relaxation
            relaxation_params =  best_candidate(result)[10:13]
        ## Kopplung
            kopplungs_params = best_candidate(result)[14:26]
        ## Degeneracy
            g = [1,2,12,16,1,3,2]#best_candidate(result)[28:end]
    # 8) Create the model

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
    Aus dem Wellenlängenscan geht hervor, dass wenn wir "rmin" pumpen, wir nicht diskret rmina pumpen, sondern auch gleichzeitig rminb. Weswegen hier ein zusätzlicher 
    Pumppulsparameter hinzugefügt wird, elcher immer 0 ist wenn rmin nicht gepumpt wird.
    Wir betrachten Modell 3.2 als Weiterentwicklung von Modell 3.1.
    Modell 3.1:
    Das Modell 3.1 ist eine Weiterentwicklung des Modells 2.2. Modell 2.2 versagt am deutlichsten bei der Beschreibung
    der rfr Mode. Hier kann keine schnelle Dynamik erreicht werden. Das Modell 3.1 versucht dies zu verbessern, indem
    es die Hin und R"uckrate zwischen der rfr Mode nicht mehr als gleich annimmt.
    Wie in Modell 2.2 beschreiben die Parameter p[1:2] die Laseranregung. Die Parameter p[3:9] beschreiben die Relaxationsraten
    der Moden. Die Parameter p[10:50] beschreiben die Kopplungsraten zwischen den Moden. Die Parameter p[51:57] beschreiben die
    Entartung der Moden.
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
    #    p[1], p[2], etc.  (Keep same indexing as your code.)
    #---------------------------------------------------------
    # Laser parameters
    α        = p[1]   # amplitude scaling of laser pulse
    α_rb     = 0   # amplitude scaling of laser pulse for rminb
    if mode == :rmin
        α_rb = p[2]
    end
    t_0      = p[3]   # center of the Gaussian pulse
    σ        = 6.67   # fixed pulse width (was hard-coded)

    # Relaxation rates (to ground state) 
    relax_rplus       = 1/p[4]
    relax_dminusomega = 1/p[5]
    relax_dfr         = 1/p[6]
    relax_dminus      = 1/p[7]
    relax_rfr         = 1/p[8]
    relax_rminusb     = 1/p[9]
    relax_rminusa     = 1/p[10]

    # Coupling constants between states (original indexing)
    # For simplicity, I'll rename them 'k_xy' to emphasize
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
    # p[31] is g_d, while g_r, g_d_omega are hard-coded to 1 below
    g_rplus   = 1 #p[38]
    g_d_omega = round(Int,p[39])
    g_dfr     = round(Int,p[40])
    g_dminus  = round(Int,p[41])
    g_rfr     = 1 #p[42]
    g_rminb   = round(Int,p[43])
    g_rmina   = round(Int,p[44])


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
    # 7) Time stepping with explicit Euler:
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
        # Control
        #-------------

        control[i,1] = control[i-1,1] + (laser_here + α_rb * f_g[i]) * dt
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )

        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        
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
    # 9) Dictionary aufbauen (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1],alpha_rb = p[2], t_0 = p[3])

    # B) Relaxation
    #    p[3..9], 7 Parameter => rplus, dminusomega, dfr, dminus, rfr, rminusb, rminusa
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[4:10]
    )

    # C) Entartung
    #    p[37..43] => 7 Parameter
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G#p[38:44]
    )

    # D) Kopplungen (jetzt mit asymmetrischen Einträgen)
    #    Wir müssen nun Vorwärts- und Rückwärts-Parameter explizit angeben.
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

    # pro Mode eine gefilterte Kopplungstabelle
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

    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Gesamtdictionary
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

    # ------------------------------------------------------
    # 10) Rückgabe als ReservoirModel (inkl. Dictionary)
    # ------------------------------------------------------
    # Falls dein ReservoirModel-Struct kein Dictionary-Feld hat,
    # kannst du hier natürlich nur das Modell und das Dict getrennt zurückgeben.
    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end

function blackbox_model_3_2(
    time_dmin,traces_dmin,use_dmin,
    time_rfr,traces_rfr,use_rfr,
    time_rmin,traces_rmin,use_rmin,
    coha;
    dt=0.1
)
    # 1) Define priors
    initial_params = [    
        5.6553,  # α_dmin
        3.4e-6,  # α_rb_dmin
        -1.98,   # t_0_dmin
        0.54,    # α_rfr
        2.52e-6, # α_rb_rfr
        -3.33,   # t_0_rfr
        6.11e-7, # α_rmin
        7.7,     # α_rb_rmin
        -2.99,   # t_0_rmin
        5000.0, # relax_rplus
        5000.0, # relax_dminusomega
        1.1594,  # relax_dfr
        1.0,     # relax_dminus
        5000.0, # relax_rfr
        5000.0, # relax_rminusb
        5000.0, # relax_rminusa
        350.68,  # rplus_dminusomega
        5000.0, # rplus_dfr
        5000.0, # rplus_dminus
        5000.0, # rplus_rfr
        35.411,  # rplus_rminusb
        5000.0, # rplus_rminusa
        5000.0, # dminusomega_dfr
        5000.0, # dminusomega_dminus
        5000.0, # dminusomega_rfr
        5000.0, # dminusomega_rminusb
        676.87,  # dminusomega_rminusa
        28.737,  # dfr_dminus
        5000.0, # dfr_rfr
        5000.0, # dfr_rminusb
        5000.0, # dfr_rminusa
        3.2224,  # dminus_rfr
        5000.0, # dminus_rminusb
        5000.0, # dminus_rminusa
        143.74,  # rfr_rplus
        7.7205,  # rfr_dminusomega
        7.6201,  # rfr_dfr
        5000.0, # rfr_dminus
        3.2201,  # rfr_rminusb
        36.081,  # rfr_rminusa
        1.0,     # rminusb_rfr
        75.208,  # rminusb_rminusa
        101.36,  # rminusa_rfr
        1.0,     # g_rplus
        1.7852,  # g_dminomega
        7.8579,  # g_dfr
        8.2365,  # g_dminus
        1.0,     # g_rfr
        2.0101,  # g_rminusb
        2.0,     # g_rminusa 
    ]
    # 2) Define the search space
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
    # 3) Construct time points for each mode
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

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
        model_dmin = model3_2(
            t_model_dmin,
            vcat(α_dmin,α_rb_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_dmin_used = model_dmin.Bleach[:, use_dmin]
        model_rfr = model3_2(
            t_model_rfr,
            vcat(α_rfr, α_rb_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        model_rmin = model3_2(
            t_model_rmin,
            vcat(α_rmin,α_rb_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rmin_used = model_rmin.Bleach[:, use_rmin]
        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
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
        # 5) Select time steps without coherence artifact in rfr
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

       
        
        # 6) Data likelihood # 6) Weights
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
        length_dmin = length(time_dmin)
        length_rfr = length(time_rfr)
        length_rmin = length(time_rmin)
        norm_weight = [1/length_dmin,1/length_rfr,1/length_rmin]
        norm_weight ./= maximum(norm_weight)
        cost_dmin = norm_weight[1] * sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
        cost_rfr =  norm_weight[2] * sum(w_rfr .* (traces_rfr - pred_rfr).^2)
        cost_rmin = norm_weight[3] * sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        # 7) Penaltys
        max_P = 0.1
        # Squared sum over all Poplulation timesteps which are above the threshold
        P_penalty = norm_weight[1] * sum(max.(model_dmin.Population .- max_P, 0.0).^2) +
            norm_weight[2] * sum(max.(model_rfr.Population .- max_P, 0.0).^2) +
            norm_weight[3] * sum(max.(model_rmin.Population .- max_P, 0.0).^2)
        total_cost = cost_dmin + cost_rfr + cost_rmin + P_penalty


        return total_cost
    end
    println("Starting optimization...")

    result = bboptimize(
        cost, # cost functionw_
        initial_params;          # initial guess
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)],
        Method       = :adaptive_de_rand_1_bin_radiuslimited ,   # or :cmaes, :genetic, etc.
        #NThreads=Threads.nthreads()-1,
        MaxSteps     = 10_000_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )

    #println("Finished. Best fitness = ", minimum_fitness(result))
    println("Best solution = ", best_candidate(result))
    # 7) Extract the parameters
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
    # 8) Create the model

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
    Das Modell 3.1 ist eine Weiterentwicklung des Modells 2.2. Modell 2.2 versagt am deutlichsten bei der Beschreibung
    der rfr Mode. Hier kann keine schnelle Dynamik erreicht werden. Das Modell 3.1 versucht dies zu verbessern, indem
    es die Hin und R"uckrate zwischen der rfr Mode nicht mehr als gleich annimmt.
    Wie in Modell 2.2 beschreiben die Parameter p[1:2] die Laseranregung. Die Parameter p[3:9] beschreiben die Relaxationsraten
    der Moden. Die Parameter p[10:50] beschreiben die Kopplungsraten zwischen den Moden. Die Parameter p[51:57] beschreiben die
    Entartung der Moden.
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
    σ        = 6.67   # fixed pulse width (was hard-coded)

    # Relaxation rates (to ground state) 
    relax_rplus       = 1/p[3]
    relax_dminusomega = 1/p[4]
    relax_dfr         = 1/p[5]
    relax_dminus      = 1/p[6]
    relax_rfr         = 1/p[7]
    relax_rminusb     = 1/p[8]
    relax_rminusa     = 1/p[9]

    # Coupling constants between states (original indexing)
    # For simplicity, I'll rename them 'k_xy' to emphasize
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
    # p[31] is g_d, while g_r, g_d_omega are hard-coded to 1 below
    g_rplus   = p[37]
    g_d_omega = p[38]
    g_dfr     = p[39]
    g_dminus  = p[40]
    g_rfr     = p[41]
    g_rminb   = p[42]
    g_rmina   = p[43]


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
    # 7) Time stepping with explicit Euler:
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
        # Control
        #-------------

        control[i,1] = control[i-1,1] + laser_here * dt
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )

        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        
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
    # 9) Dictionary aufbauen (DataFrames, etc.)
    # ------------------------------------------------------

    # A) Puls
    df_puls = DataFrame(alpha = p[1], t_0 = p[2])

    # B) Relaxation
    #    p[3..9], 7 Parameter => rplus, dminusomega, dfr, dminus, rfr, rminusb, rminusa
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode    = modes_relax,
        Wert_ps = p[3:9]
    )

    # C) Entartung
    #    p[37..43] => 7 Parameter
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = G
    )

    # D) Kopplungen (jetzt mit asymmetrischen Einträgen)
    #    Wir müssen nun Vorwärts- und Rückwärts-Parameter explizit angeben.
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

    # pro Mode eine gefilterte Kopplungstabelle
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

    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(:All => df_koppl_all)
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Gesamtdictionary
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

    # ------------------------------------------------------
    # 10) Rückgabe als ReservoirModel (inkl. Dictionary)
    # ------------------------------------------------------
    # Falls dein ReservoirModel-Struct kein Dictionary-Feld hat,
    # kannst du hier natürlich nur das Modell und das Dict getrennt zurückgeben.
    return ReservoirModel(Popu, Bleach, p, t, mode,control, outputDict)
end


function blackbox_model_3_1(
    time_dmin,traces_dmin,use_dmin,
    time_rfr,traces_rfr,use_rfr,
    time_rmin,traces_rmin,use_rmin,
    coha;
    dt=0.1
)
    # 1) Define priors
    initial_params = [    
        4.5,   # α_dmin
        -5.11,  # t_0_dmin
        0.28,   # α_rfr
        -3.03,  # t_0_rfr
        0.284,  # α_rmin
        -5.34,  # t_0_rmin
        19630,  # relax_rplus
        19809,  # relax_dminusomega
        2552,  # relax_dfr
         1.76,  # relax_dminus
         7.71,  # relax_rfr
        17933,  # relax_rminusb
        18927,  # relax_rminusa
        17.60,  # rplus_dminusomega
        19864,  # rplus_dfr
        16760,  # rplus_dminus
        19847,  # rplus_rfr
        17267,  # rplus_rminusb
        26.53,  # rplus_rminusa
        181.8,  # dminusomega_dfr
        19845,  # dminusomega_dminus
        19899,  # dminusomega_rfr
        41.81,  # dminusomega_rminusb
        14862,  # dminusomega_rminusa
        403.2,  # dfr_dminus
        19977,  # dfr_rfr
        19940,  # dfr_rminusb
        19921,  # dfr_rminusa
        31.49,  # dminus_rfr
        19981,  # dminus_rminusb
        19951,  # dminus_rminusa
        17.12,  # rfr_rplus
        19755,  # rfr_dminusomega
        19998,  # rfr_dfr
        19536,  # rfr_dminus
        19037,  # rfr_rminusb
         8.73,  # rfr_rminusa
        19133,  # rminusb_rfr
         1.00,  # rminusb_rminusa
         8.15,  # rminusa_rfr
         1.77,  # g_rplus
         7.38,  # g_dminomega
        16.98,  # g_dfr
         1.00,  # g_dminus
         1.00,  # g_rfr
         2.90,  # g_rminusb
         1.27,  # g_rminusa 
    ]
    # 2) Define the search space
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
    # 3) Construct time points for each mode
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

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

        # 2) Construct the model
        model_dmin = model3_1(
            t_model_dmin,
            vcat(α_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_dmin_used = model_dmin.Bleach[:, use_dmin]
        model_rfr = model3_1(
            t_model_rfr,
            vcat(α_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        model_rmin = model3_1(
            t_model_rmin,
            vcat(α_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rmin_used = model_rmin.Bleach[:, use_rmin]
        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
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
        # 5) Select time steps without coherence artifact in rfr
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = time_rfr[selected_idx_rfr]

        # 6) Weights
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
        
        
        # 6) Data likelihood
        cost_dmin = sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
        cost_rfr =  sum(w_rfr .* (traces_rfr - pred_rfr).^2)
        cost_rmin = sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        total_cost = cost_dmin + cost_rfr + cost_rmin

        return total_cost
    end
    println("Starting optimization...")

    result = bboptimize(
        cost, # cost functionw_
        initial_params;          # initial guess
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)],
        Method       = :adaptive_de_rand_1_bin_radiuslimited ,   # or :cmaes, :genetic, etc.
        #NThreads=Threads.nthreads()-1,
        MaxSteps     = 10_000_000, # or some large budget
        PopulationSize = 2_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )

    #println("Finished. Best fitness = ", minimum_fitness(result))
    println("Best solution = ", best_candidate(result))
    # 7) Extract the parameters
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
    # 8) Create the model

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
    Modell 2.2:
    Beschreibt die Dynamik aller Moden, welche eine kontinuierliche "Anderung in Experimenten (ODT,UDT,HT)
    aufweisen. Die Moden sind: rplus, dminusomega, dfr, dminus, rfr, rminusb, rminusa.
    Jede Mode kann durch einen Laser mit den Parametern α und t₀ angeregt werden. Daf"ur kwarg `mode=:dmin`etc. verwenden.
    Jede Mode kann relaxieren p[3:9] oder mit den anderen Moden wechselwirken p[10:30].
    Weiterhin kann die Entartung der Moden mit p[31:37] eingestellt werden.
"""
function model2_2(
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
    σ        = 6.67   # fixed pulse width (was hard-coded)

    # Relaxation rates (to ground state) 
    relax_rplus       = 1/p[3]
    relax_dminusomega = 1/p[4]
    relax_dfr         = 1/p[5]
    relax_dminus      = 1/p[6]
    relax_rfr         = 1/p[7]
    relax_rminusb     = 1/p[8]
    relax_rminusa     = 1/p[9]

    # Coupling constants between states (original indexing)
    # For simplicity, I'll rename them 'k_xy' to emphasize
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

    k_rfr_rplus           = k_rplus_rfr
    k_rfr_dminusomega     = k_dminusomega_rfr
    k_rfr_dfr             = k_dfr_rfr
    k_rfr_dminus          = k_dminus_rfr
    k_rfr_rminusb         = 1/p[28]
    k_rfr_rminusa         = 1/p[29]

    k_rminusb_rplus       = k_rplus_rminusb
    k_rminusb_dminusomega = k_dminusomega_rminusb
    k_rminusb_dfr         = k_dfr_rminusb
    k_rminusb_dminus      = k_dminus_rminusb
    k_rminusb_rfr         = k_rfr_rminusb
    k_rminusb_rminusa     = 1/p[30]

    k_rminusa_rplus       = k_rplus_rminusa
    k_rminusa_dminusomega = k_dminusomega_rminusa
    k_rminusa_dfr         = k_dfr_rminusa
    k_rminusa_dminus      = k_dminus_rminusa
    k_rminusa_rfr         = k_rfr_rminusa
    k_rminusa_rminusb     = k_rminusb_rminusa

    # 3) Degeneracy of the states:
    # p[31] is g_d, while g_r, g_d_omega are hard-coded to 1 below
    g_rplus   = p[31]
    g_d_omega = p[32]
    g_dfr     = p[33]
    g_dminus  = p[34]
    g_rfr     = p[35]
    g_rminb   = p[36]
    g_rmina   = p[37]


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
    # 7) Time stepping with explicit Euler:
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
        # Control
        #-------------

        control[i,1] = control[i-1,1] + laser_here * dt
        control[i,2] = G[1] * Popu[i,1] + (
               G[2] * Popu[i,2] 
             + G[3] * Popu[i,3] 
             + G[4] * Popu[i,4] 
             + G[5] * Popu[i,5] 
             + G[6] * Popu[i,6] 
             + G[7] * Popu[i,7] 
        )

        control[i,3] = control[i-1,3] + (
            + G[1] * relax_rplus * Popu[i-1,1]
            + G[2] * relax_dminusomega * Popu[i-1,2]
            + G[3] * relax_dfr * Popu[i-1,3]
            + G[4] * relax_dminus * Popu[i-1,4]
            + G[5] * relax_rfr * Popu[i-1,5]
            + G[6] * relax_rminusb * Popu[i-1,6]
            + G[7] * relax_rminusa * Popu[i-1,7]
        ) * dt
        
        control[i,4] = control[i,1] - control[i,2] - control[i,3]
    end

    #---------------------------------------------------------
    # 8) Finally, compute the "Bleach" signal 
    #    for each state: (1 - 2*g[i]*Popu[i])^2
    #---------------------------------------------------------
    for k in 1:n
        Bleach[:,k] = (1 .- 2 .* G[k] .* Popu[:,k]) .^2
    end

    #---------------------------------------------------------
    # 9) Return the results in a struct or dictionary
    #---------------------------------------------------------
    # -- A) Puls-Parameter
     df_puls = DataFrame(
        alpha = p[1],
        t_0   = p[2]
    )

    # -- B) Relaxations-Parameter als DataFrame
    #    Zuordnung: p[3] -> rplus, p[4] -> dminusomega, ...
    modes_relax = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_relax = DataFrame(
        Mode = modes_relax,
        Wert_ps = p[3:9]  # p[3] bis p[9]
    )

    # -- C) Entartung (Degeneracies) als DataFrame
    #    Zuordnung: p[31] -> rplus, p[32] -> dminusomega, ...
    modes_deg = ["rplus", "dminusomega", "dfr", "dminus", "rfr", "rminusb", "rminusa"]
    df_entartung = DataFrame(
        Mode = modes_deg,
        g    = p[31:37]  # p[31] bis p[37]
    )

    # -- D) Kopplungen:
    #    Wir listen alle k_x_y auf, mit Index = [10..30].
    #    k_rplus_dminusomega = 1/p[10], also Parameter p[10], usw.
    #    Für jede Kopplung (x,y,pIndex) legen wir zwei Zeilen an:
    #       t_x_y   => p[pIndex],
    #       t_y_x   => p[pIndex].
    #    (Da rückwärts = vorwärts in deinem Modell.)
    #
    #    Wichtig: "k_" wird durch "t_" ersetzt.

    # Hilfstabelle, die eindeutig zuordnet:
    # ( :rplus, :dminusomega, 10 ) bedeutet,
    #    k_rplus_dminusomega => 1/p[10] in deinem Modell
    #    Aber wir wollen "t_rplus_dminusomega" => p[10] in DataFrame
    #    analog t_dminusomega_rplus => p[10].
    couplings = [
        (:rplus,       :dminusomega, 10),
        (:rplus,       :dfr,         11),
        (:rplus,       :dminus,      12),
        (:rplus,       :rfr,         13),
        (:rplus,       :rminusb,     14),
        (:rplus,       :rminusa,     15),
        (:dminusomega, :dfr,         16),
        (:dminusomega, :dminus,      17),
        (:dminusomega, :rfr,         18),
        (:dminusomega, :rminusb,     19),
        (:dminusomega, :rminusa,     20),
        (:dfr,         :dminus,      21),
        (:dfr,         :rfr,         22),
        (:dfr,         :rminusb,     23),
        (:dfr,         :rminusa,     24),
        (:dminus,      :rfr,         25),
        (:dminus,      :rminusb,     26),
        (:dminus,      :rminusa,     27),
        (:rfr,         :rminusb,     28),
        (:rfr,         :rminusa,     29),
        (:rminusb,     :rminusa,     30),
    ]

    # Alle Kopplungen in einer "All"-Tabelle:
    rows_all = Vector{Dict{Symbol,Any}}()
    for (x,y, idx) in couplings
        # Vorwärts
        push!(rows_all, Dict(
            :param_name => "t_$(x)_$(y)",  # z.B. "t_rplus_dminusomega"
            :value      => p[idx]
        ))
        # Rückwärts
        push!(rows_all, Dict(
            :param_name => "t_$(y)_$(x)",
            :value      => p[idx]
        ))
    end
    df_koppl_all = DataFrame(rows_all)

    # Nun zusätzlich pro Mode eine gefilterte Tabelle,
    # in der nur Kopplungen auftauchen, in denen die gegebene Mode
    # entweder x oder y ist.
    function df_koppl_for_mode(m::Symbol)
        rows = Vector{Dict{Symbol,Any}}()
        for (x,y, idx) in couplings
            if x == m || y == m
                # vorwärts
                push!(rows, Dict(
                    :param_name => "t_$(x)_$(y)",
                    :value      => p[idx]
                ))
                # rückwärts
                push!(rows, Dict(
                    :param_name => "t_$(y)_$(x)",
                    :value      => p[idx]
                ))
            end
        end
        return DataFrame(rows)
    end

    # Jetzt bauen wir ein Dict daraus, wo Key=Symbol der Mode ist und Value = DataFrame
    modes_all = [:rplus, :dminusomega, :dfr, :dminus, :rfr, :rminusb, :rminusa]
    dict_koppl = Dict(
        :All => df_koppl_all
    )
    # jede Mode einzeln ergänzen
    for m in modes_all
        dict_koppl[m] = df_koppl_for_mode(m)
    end

    # Gesamtstruktur
    outputDict = Dict(
        :fitted_mode => mode,
        :Parameter   => Dict(
            :All         => p,
            :Puls        => df_puls,
            :Relaxation  => df_relax,
            :Entartung   => df_entartung,
            :Kopplung    => dict_koppl
        ),
        :Model       => "Model 2.2",
        :Zeitschritt => dt
    )

    # Return 
    return ReservoirModel(Popu, Bleach, p, t, mode, control, outputDict)
end


function blackbox_model_2_2(
    time_dmin,traces_dmin,use_dmin,
    time_rfr,traces_rfr,use_rfr,
    time_rmin,traces_rmin,use_rmin,
    coha;
    dt=0.1
)
    # 1) Define priors
    initial_params = vcat(
        [0.13, 0.0], # α_dmin, t_0_dmin
        [0.075, 5.0], # α_rfr, t_0_rfr
        [0.2, -7.0], # α_rmin, t_0_rmin
        fill(500,7), # 7 relaxation_params   
        fill(500,21), # 21 kopplungs_params
        fill(2.0,7)
    )
    # 2) Define the search space
    lb = vcat(
        [0.0, -20.0], # α_dmin, t_0_dmin
        [0.0, -20.0], # α_rfr, t_0_rfr
        [0.0, -20.0], # α_rmin, t_0_rmin
        fill(1.0,7), # relaxation_params   
        fill(1.0,21), # kopplungs_params
        fill(1.0,7) # g_params
    )
    ub = vcat(
        [1.0, 20.0], # α_dmin, t_0_dmin
        [3.0, 20.0], # α_rfr, t_0_rfr
        [1.0, 20.0], # α_rmin, t_0_rmin
        fill(20000.0,7), # relaxation_params   
        fill(20000.0,21), # kopplungs_params
        fill(17.0,7) # g_params
    )
    # 3) Construct time points for each mode
    t_model_dmin = collect(time_dmin[1]:dt:time_dmin[end])
    t_model_rfr = collect(time_rfr[1]:dt:time_rfr[end])
    t_model_rmin = collect(time_rmin[1]:dt:time_rmin[end])

    function cost(
        p::Vector{Float64}
    )
        # 1) Unpack the parameter vector
        α_dmin  = p[1];  t0_dmin  = p[2]
        α_rfr   = p[3];  t0_rfr   = p[4]
        α_rmin  = p[5];  t0_rmin  = p[6]

        relax_p = p[7:13]      # 7 relaxation
        koppl_p = p[14:34]     # 21 couplings
        g     = p[35:end]      # 5 degeneracy

        # 2) Construct the model
        model_dmin = model2_2(
            t_model_dmin,
            vcat(α_dmin, t0_dmin, relax_p, koppl_p, g),
            mode = :dmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_dmin_used = model_dmin.Bleach[:, use_dmin]
        model_rfr = model2_2(
            t_model_rfr,
            vcat(α_rfr, t0_rfr, relax_p, koppl_p, g),
            mode = :rfr,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rfr_used = model_rfr.Bleach[:, use_rfr]
        model_rmin = model2_2(
            t_model_rmin,
            vcat(α_rmin, t0_rmin, relax_p, koppl_p, g),
            mode = :rmin,
            dt=dt
        )
        # We'll use only a subset of the excited states since not all modes show a continous change
        model_rmin_used = model_rmin.Bleach[:, use_rmin]
        # 3) Interpolate the predicted signals at the experiment time points
        spl_dmin = [linear_interpolation(t_model_dmin, model_dmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_dmin_used,2)]
        spl_rfr  = [linear_interpolation(t_model_rfr, model_rfr_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rfr_used,2)]
        spl_rmin = [linear_interpolation(t_model_rmin, model_rmin_used[:,i], extrapolation_bc=Line()) for i in 1:size(model_rmin_used,2)]

        v_pred_dmin =  [spl_bleach(time_dmin) for spl_bleach in spl_dmin]
        v_pred_rfr  =  [spl_bleach(time_rfr) for spl_bleach in spl_rfr]
        v_pred_rmin =  [spl_bleach(time_rmin) for spl_bleach in spl_rmin]

        # 4) Combine them to match data's shape.
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
        # 5) Select time steps without coherence artifact in rfr
        selected_idx_rfr = setdiff(1:size(time_rfr,1),last(coha)[1])
        selected_time_rfr = deepcopy(time_rfr[selected_idx_rfr])
        # 6) Weights
        w_dmin = ones(size(pred_dmin))
        w_pertrace_dmin = [
            5.0,  # rplus
            1.0,  # rfr
            1.0,  # rminusb
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
            1.0,  # rfr
            1.0,  # rminusb
            5.0   # rminusa
        ]
        for m in 1:size(w_rmin,1)
            for n in 1:size(w_rmin,2)
                w_rmin[m,n] = time_weight(time_rmin[m],weight=3) * w_pertrace_rmin[n]
            end
        end
        
        # 6) Data likelihood
        cost_dmin = sum(w_dmin .* (traces_dmin - pred_dmin).^2 )
        cost_rfr =  sum(w_rfr  .* (traces_rfr - pred_rfr).^2)
        cost_rmin = sum(w_rmin .* (traces_rmin - pred_rmin).^2)
        total_cost = cost_dmin + cost_rfr + cost_rmin

        return total_cost
    end
    println("Starting optimization...")

    result = bboptimize(
        cost, # cost functionw_
        #initial_params;          # initial guess
        [0.42931043246275824, -8.146, 0.2863252315666348, -6.6, 0.32853329837264816, -5.572319999012165, 845.1397263180909, 181.77299824914687, 364.47593379263895, 51.88005710270164, 1731.8185341148842, 692.6863710034145, 382.27783620987833, 1260.8989693159294, 423.07198576473127, 1188.2549219347056, 2066.7479200695866, 1.0898267697002924, 1358.1301625942042, 1.04165287690823966, 1887.2648541668257, 720.0502887648229, 1217.8148634149704, 743.2048836079019, 951.4914851184341, 14.89227367991587, 1582.2929493340866, 224.00668662027235, 49.359739628032216, 311.5398513802186, 180.13303255017343, 1.00990611155797902, 1.00110428169553246, 149.12608088781172, 2.3867006900820544, 11.335361158399444, 16.99947051974148, 5.396202493170724, 1.0000003220906513, 4.143117000442034, 1.7180100885473548],
        #[0.451511, 0.155268, 0.320754, 8.95625, 0.464752, -3.83644, 34.7093, 1235.34, 397.267, 940.581, 2086.29, 836.934, 27.9323, 0.105189, 1460.91, 601.86, 853.731, 775.261, 567.368, 1479.01, 2155.82, 1817.02, 0.103642, 874.106, 0.16992, 927.949, 26.9699, 273.941, 1590.06, 1248.41, 2351.63, 0.106644, 0.361485, 1044.54, 2.09593, 7.48442, 17.0, 16.9882, 1.00408, 3.21444, 1.44737],
        #[0.402576, -15.0, 0.306215, 14.1604, 0.391271, -4.90299, 350.022, 119.957, 640.174, 267.288, 732.974, 1491.4, 549.259, 1.93351, 1.25758, 142.866, 660.834, 1006.69, 593.691, 762.518, 577.137, 896.104, 938.247, 1469.86, 487.144, 5.15297, 630.948, 266.348, 420.879, 487.044, 733.371, 1.41609, 690.943, 1.0, 2.34514, 9.05271, 17.0, 17.0, 1.0, 3.64803, 1.48169],
        SearchRange = [(lb[i], ub[i]) for i in 1:length(lb)],
        Method       = :dxnes  ,   # or :cmaes, :genetic, etc.
        #NThreads=Threads.nthreads()-1,
        MaxSteps     = 50_000, # or some large budget
        PopulationSize = 1_000,  # can try bigger for better coverage
        TraceMode    = :verbose 
    )

    #println("Finished. Best fitness = ", minimum_fitness(result))
    println("Best solution = ", best_candidate(result))
    # 7) Extract the parameters
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
            kopplungs_params = best_candidate(result)[14:34]
        ## Degeneracy
            g = best_candidate(result)[35:end]
    # 8) Create the model

    best_model_dmin  = model2_2(
        t_model_dmin,
        vcat(pulse_params_dmin,relaxation_params,kopplungs_params,g),
        mode=:dmin,
        dt=dt
    )
    best_model_rfr  = model2_2(
        t_model_rfr,
        vcat(pulse_params_rfr,relaxation_params,kopplungs_params,g),
        mode=:rfr,
        dt=dt
    )
    best_model_rmin  = model2_2(
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
"""
function time_weight(t;
    weight = 3.0
    )
    if -15.0 <= t < 30.0
        return weight  # heavier weight for early times
    else
        return 1.0   # normal weight for later times
    end
end