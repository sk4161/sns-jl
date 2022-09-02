using PyPlot

const TAU = 20.0                    # ms
const V_REST = -65.0                # mV
const V_RESET = -65.0               # mV
const THETA = -55.0                 # mV
const R_M = 1.0                     # MΩ
const DT = 1.0                      # ms
const T = 1000                      # ms
const NT = T / DT        
const I_EXT = 12.0                  # nA
const N = 2                         # of neurons

function euler()
    v = [V_REST, V_REST - 15.0]
    s = [false, false]

    t_list = [[], []]
    v_list = [[], []]

    for nt = 0:NT-1  
        t = DT * nt

        for i = 1:N
            push!(t_list[i], t)
            push!(v_list[i], v[i])
        end

        for i = 1:N
            v[i] += DT * (-(v[i] - V_REST) + R_M * I_EXT) / TAU
            s[i] = (v[i] > THETA)
        end

        # Pretty-print spikes on membrane potentials. Note that spike time is not t but t + DT
        
        for i = 1:N
            if s[i]
                push!(t_list[i], t + DT)        
                push!(v_list[i], v[i])
                push!(t_list[i], t + DT)
                push!(v_list[i], 0.0)
            end
        end

        for i = 1:N
            v[i] = s[i] * V_RESET + (!s[i]) * v[i]
        end

    end

    title("LIF model (2 neurons)")
    xlabel("Times [ms]")
    ylabel("Membrane Potential [mV]")
    xlim(0, T)
    ylim(-70, 0)
    for i = 1:N
        plot(t_list[i], v_list[i], label="Neuron $i")
    end
    legend()
    savefig("../../../fig/part1/lif/lif2.png")
end

euler()