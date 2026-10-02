"""
    REMPILookup

Optional detection model for `IntDist.product_spectrum`, in which the ionised and
surviving fractions come from integrating a three-level ionisation ODE over the REMPI
pulse rather than from a linear intensity profile. Kept separate because it needs
DifferentialEquations, which is slow to load.

Include after IntDist.jl:

    include("REMPILookup.jl"); using .REMPILookup
    det = lookup_detection(rempi_cfg)
"""
module REMPILookup

using DifferentialEquations
using Interpolations
import ..IntDist: REMPIConfig

export rempi_lookup, lookup_detection

# y = (ground amplitude, intermediate amplitude, accumulated ion population).
# Standard Rabi-ionisation form: dc_g = i(Ω/2)c_e, dc_e = i(Ω/2)c_g - (Γ/2)c_e,
# dP_ion = Γion|c_e|^2. The original notebook had dy[1] = -0.5im*(Ω*y[2] - conj(Ω*y[2])),
# which does not conserve population (survival up to ~22, ion yield up to ~3.6).
function rempi_ode!(dy, y, p, t)
    Ω, Γion, Γdecay = p
    dy[1] = 0.5im * conj(Ω) * y[2]
    dy[2] = -0.5 * (Γdecay + Γion) * y[2] + 0.5im * Ω * y[1]
    dy[3] = Γion * abs2(y[2])
end

"""
    rempi_lookup(max_rabi, pulse_width, Γion, Γdecay; n = 100)

Tabulate the ion yield and ground-state survival after one pulse against Rabi
frequency (rad/s) from 0 to `max_rabi`, and return linear interpolants of both.
"""
function rempi_lookup(max_rabi, pulse_width, Γion, Γdecay; n = 100)
    omegas = range(0, max_rabi, length = n)
    ion  = zeros(n)
    surv = ones(n)
    for (i, Ω) in enumerate(omegas)
        Ω == 0 && continue
        prob = ODEProblem(rempi_ode!, ComplexF64[1, 0, 0], (0.0, pulse_width), (Ω, Γion, Γdecay))
        y = DifferentialEquations.solve(prob, Tsit5(), save_everystep = false, verbose = false).u[end]
        ion[i]  = real(y[3])
        surv[i] = abs2(y[1])
    end
    return linear_interpolation(omegas, ion, extrapolation_bc = Flat()),
           linear_interpolation(omegas, surv, extrapolation_bc = Flat())
end

"""
    lookup_detection(cfg; max_rabi = 2π*15e6, pulse_width = 50e-9,
                     Γion = 2π*14e6, Γdecay = 2π*14e6)

Detection model `r -> (ionised, surviving)` for `product_spectrum`. Outside the dark
line the local Rabi frequency is `max_rabi * exp(-r^2 / w^2)` (field, not intensity,
profile); inside it is zero. `cfg.max_ionization_prob` is not used.
"""
function lookup_detection(cfg::REMPIConfig; max_rabi = 2pi * 15e6, pulse_width = 50e-9,
                          Γion = 2pi * 14e6, Γdecay = 2pi * 14e6)
    itp_ion, itp_surv = rempi_lookup(max_rabi, pulse_width, Γion, Γdecay)
    return function (r)
        Ω = r >= cfg.dark_line_radius ? max_rabi * exp(-r^2 / cfg.beam_waist^2) : 0.0
        return itp_ion(Ω), itp_surv(Ω)
    end
end

end # module
