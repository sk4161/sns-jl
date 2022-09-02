using PyPlot

const E_REST = -65.0            # mV
const C = 1.0                   # μF / cm^2
const G_LEAK = 0.3              # mS / cm^2
const E_LEAK = 10.6 + E_REST    # mV 
const G_NA = 120.0              # mS / cm^2
const E_NA = 115.0 + E_REST     # mV
const G_K = 36.0                # mS / cm^2
const E_K = -12.0 + E_REST      # mV

const DT = 0.01                 # ms
const T = parse(Int64, ARGS[1])
const NT = T / DT               # 1000 ms / DT

alpha_m(v) = (2.5 - 0.1 * (v - E_REST)) / (exp.(2.5 - 0.1 * (v - E_REST)) - 1.0)
beta_m(v) = 4.0 * exp.(-(v - E_REST) / 18.0)
alpha_h(v) = 0.07 * exp.(-(v - E_REST) / 20.0)
beta_h(v) = 1.0 / (exp.(3.0 - 0.1 * (v - E_REST)) + 1.0)
alpha_n(v) = (0.1 - 0.01 * (v - E_REST)) / (exp.(1.0 - 0.1 * (v - E_REST)) - 1.0)
beta_n(v) = 0.125 * exp.(-(v - E_REST) / 80.0)

m0(v) = alpha_m(v) / (alpha_m(v) + beta_m(v))
h0(v) = alpha_h(v) / (alpha_h(v) + beta_h(v))
n0(v) = alpha_n(v) / (alpha_n(v) + beta_n(v))
tau_m(v) = 1.0 / (alpha_m(v) + beta_m(v))
tau_h(v) = 1.0 / (alpha_h(v) + beta_h(v))
tau_n(v) = 1.0 / (alpha_n(v) + beta_n(v))

dmdt(v, m) = (1.0 / tau_m(v)) * (-m + m0(v))
dhdt(v, h) = (1.0 / tau_h(v)) * (-h + h0(v))
dndt(v, n) = (1.0 / tau_n(v)) * (-n + n0(v))
dvdt(v, m, h, n, i_ext) = (-G_LEAK * (v - E_LEAK) - G_NA * m^3 * h * (v - E_NA) - G_K * n^4 * (v - E_K) + i_ext) / C

function runge_kutta()
    v = E_REST
    m = m0(v)
    h = h0(v)
    n = n0(v)

    i_ext = 9.0 # μA / cm^2

    t_list = []
    v_list = []

    for nt = 0:NT-1
        t = DT * nt

        push!(t_list, t)
        push!(v_list, v)

        dmdt1 = dmdt(v, m)
        dhdt1 = dhdt(v, h)
        dndt1 = dndt(v, n)
        dvdt1 = dvdt(v, m, h, n, i_ext)

        dmdt2 = dmdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1)
        dhdt2 = dhdt(v + 0.5 * DT * dvdt1, h + 0.5 * DT * dhdt1)
        dndt2 = dndt(v + 0.5 * DT * dvdt1, n + 0.5 * DT * dndt1)
        dvdt2 = dvdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1, h + 0.5 * DT * dhdt1, n + 0.5 * DT * dndt1, i_ext)

        dmdt3 = dmdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2)
        dhdt3 = dhdt(v + 0.5 * DT * dvdt2, h + 0.5 * DT * dhdt2)
        dndt3 = dndt(v + 0.5 * DT * dvdt2, n + 0.5 * DT * dndt2)
        dvdt3 = dvdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2, h + 0.5 * DT * dhdt2, n + 0.5 * DT * dndt2, i_ext)

        dmdt4 = dmdt(v + DT * dvdt3, m + DT * dmdt3)
        dhdt4 = dhdt(v + DT * dvdt3, h + DT * dhdt3)
        dndt4 = dndt(v + DT * dvdt3, n + DT * dndt3)
        dvdt4 = dvdt(v + DT * dvdt3, m + DT * dmdt3, h + DT * dhdt3, n + DT * dndt3, i_ext)

        m += DT * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6.0
        h += DT * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6.0
        n += DT * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6.0
        v += DT * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6.0
    end

    suptitle("HH model")
    subplot(211, title="$(T/10) ms", xlim=(0, T/10), ylim=(-80, 60), xlabel="Times [ms]", ylabel="Membrane Potential [mV]")
    plot(t_list, v_list)
    subplot(212, title="$(T) ms", xlim=(0, T), ylim=(-80, 60), xlabel="Times [ms]", ylabel="Membrane Potential [mV]")
    plot(t_list, v_list)
    tight_layout()

    savefig("../../../fig/part1/hh/hh.png")
end

runge_kutta()