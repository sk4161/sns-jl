using PyPlot

const TAU = 20.0                    # ms
const V_REST = -65.0                # mV
const V_RESET = -65.0               # mV
const THETA = -55.0                 # mV
const R_M = 1.0                     # MΩ
const DT = 1.0                      # ms
const T = parse(Int64, ARGS[1])     # ms
const NT = T / DT        
const I_EXT = 12.0                  # nA
const T_REFR = 5.0                  # ms
const NT_REFR = T_REFR / DT

function euler()
    v = V_REST
    refr = 0

    t_list = []
    v_list = []

    for nt = 0:NT-1  
        t = DT * nt

        push!(t_list, t)
        push!(v_list, v)

        v += DT * (-(v - V_REST) + R_M * I_EXT) / TAU
        s = (v > THETA)

        # Pretty-print spikes on membrane potentials. Note that spike time is not t but t + DT
        if s        
            push!(t_list, t + DT)
            push!(v_list, v)
            push!(t_list, t + DT)
            push!(v_list, 0.0)
        end
        refr = s * NT_REFR + (!s) * (refr - 1)
        v = (refr > 0) ? V_RESET : v
    end

    title("LIF model ($T ms)")
    xlabel("Times [ms]")
    ylabel("Membrane Potential [mV]")
    xlim(0, T)
    ylim(-70, 0)
    plot(t_list, v_list)
    savefig("../../../fig/part1/lif/lif_refr_$(T)_ms.png")
end

euler()