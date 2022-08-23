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

function euler()
    v = V_REST

    t_list = []
    v_list = []

    for nt = 0:NT-1  
        t = DT * nt

        push!(t_list, t)
        push!(v_list, v)

        s = (v > THETA)
        dv = DT * (-(v - V_REST) + R_M * I_EXT) / TAU

        v = s * V_RESET + (!s) * (v + dv)
    end

    title("LIF model ($T ms)")
    xlabel("Times [ms]")
    ylabel("Membrane Potential [mV]")
    xlim(0, T)
    ylim(-65, -55)
    plot(t_list, v_list)
    savefig("../../../fig/part1/lif/lif_alt_$(T)_ms.png")
end

euler()