using PyPlot

const C = 1.0                   # μF / cm^2
const G_LEAK = 0.3              # mS / cm^2
const E_LEAK = -17.0            # mV 
const G_NA = 120.0              # mS / cm^2
const E_NA = 55.0               # mV
const G_K = 20.0                # mS / cm^2
const E_K = -75.0               # mV
const G_A = 47.7                # mS / cm^2
const E_A = -75.0               # mV

const DT = 0.01                 # ms
const T = 1000
const NT = T / DT               # 1000 ms / DT

const I_EXT = parse(Float64, ARGS[1]) # μA / cm^2

const mshft = -5.3
const hshft = -12.0
const nshft = -4.3

alpha_m(v) = (-0.1 * (v + 35.0 + mshft)) / (exp.(- 0.1 * (v + 35.0 + mshft)) - 1.0)
beta_m(v) = 4.0 * exp.(-(v + 60.0 + mshft) / 18.0)
alpha_h(v) = 0.07 * exp.(-(v + 60.0 + hshft) / 20.0)
beta_h(v) = 1.0 / (exp.(-0.1 * (v + 30.0 + hshft)) + 1.0)
alpha_n(v) = (- 0.01 * (v + 50.0 + nshft)) / (exp.(- 0.1 * (v + 50.0 + nshft)) - 1.0)
beta_n(v) = 0.125 * exp.(-(v + 60 + nshft) / 80.0)

m0(v) = alpha_m(v) / (alpha_m(v) + beta_m(v))
h0(v) = alpha_h(v) / (alpha_h(v) + beta_h(v))
n0(v) = alpha_n(v) / (alpha_n(v) + beta_n(v))
a0(v) = cbrt.(0.0761 * exp.((v + 94.22) / 31.84) / (1.0 + exp.((v + 1.17) / 28.93)))
b0(v) = 1.0 / (1.0 + exp.((v + 53.3) / 14.54))^4
tau_m(v) = (1.0 / 3.8) * (1.0 / (alpha_m(v) + beta_m(v)))
tau_h(v) = (1.0 / 3.8) * (1.0 / (alpha_h(v) + beta_h(v)))
tau_n(v) = 2.0 * (1.0 / (alpha_n(v) + beta_n(v)))
tau_a(v) = 0.3632 + 1.158 / (1.0 + exp.((v + 55.96) / 20.12))
tau_b(v) = 1.24 + 2.678 / (1.0 + exp.((v + 50.0) / 16.027))

dmdt(v, m) = (1.0 / tau_m(v)) * (-m + m0(v))
dhdt(v, h) = (1.0 / tau_h(v)) * (-h + h0(v))
dndt(v, n) = (1.0 / tau_n(v)) * (-n + n0(v))
dadt(v, a) = (1.0 / tau_a(v)) * (-a + a0(v))
dbdt(v, b) = (1.0 / tau_b(v)) * (-b + b0(v))
dvdt(v, m, h, n, a, b, i_ext) = (- G_LEAK * (v - E_LEAK) - G_NA * m^3 * h * (v - E_NA) - G_K * n^4 * (v - E_K)
                                 - G_A * a^3 * b * (v - E_A) # I_A-current
                                 + i_ext) / C

function runge_kutta()
    v = E_LEAK
    m = m0(v)
    h = h0(v)
    n = n0(v)
    a = a0(v)
    b = b0(v)

    t_list = []
    v_list = []

    for nt = 0:NT-1
        t = DT * nt

        push!(t_list, t)
        push!(v_list, v)

        dmdt1 = dmdt(v, m)
        dhdt1 = dhdt(v, h)
        dndt1 = dndt(v, n)
        dadt1 = dadt(v, a)
        dbdt1 = dbdt(v, b)
        dvdt1 = dvdt(v, m, h, n, a, b, I_EXT)

        dmdt2 = dmdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1)
        dhdt2 = dhdt(v + 0.5 * DT * dvdt1, h + 0.5 * DT * dhdt1)
        dndt2 = dndt(v + 0.5 * DT * dvdt1, n + 0.5 * DT * dndt1)
        dadt2 = dadt(v + 0.5 * DT * dvdt1, a + 0.5 * DT * dadt1)
        dbdt2 = dbdt(v + 0.5 * DT * dvdt1, b + 0.5 * DT * dbdt1)
        dvdt2 = dvdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1, h + 0.5 * DT * dhdt1, n + 0.5 * DT * dndt1, 
                     a + 0.5 * DT * dadt1, b + 0.5 * DT * dbdt1, I_EXT)

        dmdt3 = dmdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2)
        dhdt3 = dhdt(v + 0.5 * DT * dvdt2, h + 0.5 * DT * dhdt2)
        dndt3 = dndt(v + 0.5 * DT * dvdt2, n + 0.5 * DT * dndt2)
        dadt3 = dadt(v + 0.5 * DT * dvdt2, a + 0.5 * DT * dadt2)
        dbdt3 = dbdt(v + 0.5 * DT * dvdt2, b + 0.5 * DT * dbdt2)
        dvdt3 = dvdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2, h + 0.5 * DT * dhdt2, n + 0.5 * DT * dndt2, 
                     a + 0.5 * DT * dadt2, b + 0.5 * DT * dbdt2, I_EXT)

        dmdt4 = dmdt(v + DT * dvdt3, m + DT * dmdt3)
        dhdt4 = dhdt(v + DT * dvdt3, h + DT * dhdt3)
        dndt4 = dndt(v + DT * dvdt3, n + DT * dndt3)
        dadt4 = dadt(v + DT * dvdt3, a + DT * dadt3)
        dbdt4 = dbdt(v + DT * dvdt3, b + DT * dbdt3)
        dvdt4 = dvdt(v + DT * dvdt3, m + DT * dmdt3, h + DT * dhdt3, n + DT * dndt3, 
                     a + DT * dadt3, b + DT * dbdt3, I_EXT)

        m += DT * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6.0
        h += DT * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6.0
        n += DT * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6.0
        a += DT * (dadt1 + 2 * dadt2 + 2 * dadt3 + dadt4) / 6.0
        b += DT * (dbdt1 + 2 * dbdt2 + 2 * dbdt3 + dbdt4) / 6.0
        v += DT * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6.0
    end

    title("Connor-Stevens model ($(I_EXT) " * raw"$\mu$A/$\mathrm{cm^2})$")
    xlabel("Times [ms]")
    ylabel("Membrane Potential [mV]")
    xlim(0, T)
    ylim(-80, 60)
    plot(t_list, v_list)
    savefig("../../../fig/part1/hh/ia_$(I_EXT).png")
end

runge_kutta()