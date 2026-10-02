"""
    IntDist

Microwave spectroscopy of KRb hyperfine states in an optical dipole trap.

Given a set of populated N=0 hyperfine states, simulate the spectrum seen when a
microwave field is scanned across the N=0 -> N=1 transitions, either at a single
trap intensity or averaged over the thermal intensity distribution sampled by the
molecules (the "IntDist" in the name).

Pipeline (see IntDist.ipynb):

 1. Trap geometry and thermal sampling -> per-molecule intensities
      `GaussianBeam`, `trap_intensity`, `trap_potential`, `sample_thermal`
 2. Molecular structure at a grid of intensities -> `DressedMolecule`
      (built from a `diatomic_jl` `solve.sol` plus the dipole operators)
 3. Microwave drive -> `MicrowaveDrive`
 4. Spectra
      `spectrum`                     full N<=1 space, one intensity
      `spectrum_truncated`           energy-window truncation, one intensity
      `intensity_averaged_spectrum`  truncated spectrum averaged over molecules
      `multistep_spectrum`           repeated pulse + partial measurement
      `product_spectrum`             Monte Carlo REMPI detection of flying products

Energies are in Hz throughout (not angular), times in s, intensities in W/m^2.
"""
module IntDist

using LinearAlgebra
using Random
using Statistics
using Base.Threads

export kB, h, amu, Debye, KRb_mass, KRb_polarizability_1064,
       GaussianBeam, beam_profile, trap_intensity, trap_potential,
       sample_thermal, thermal_stats,
       spherical_components, rotation_matrix, MicrowaveDrive,
       manifolds, DressedMolecule, rotating_frame_hamiltonian,
       spectrum, spectrum_truncated, intensity_averaged_spectrum,
       multistep_spectrum,
       product_speeds, REMPIConfig, linear_detection, product_spectrum

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

const kB    = 1.3806452e-23
const h     = 6.62607015e-34
const amu   = 1.66054e-27
const Debye = 3.33564e-30

const KRb_mass = (40 + 87) * amu

"KRb ground-state polarizability at 1064 nm, 55.3 Hz/(W/cm^2), in J/(W/m^2)."
const KRb_polarizability_1064 = 5.53e-5 * 1e6 * h * 1e-4

# ---------------------------------------------------------------------------
# Trap geometry
# ---------------------------------------------------------------------------

"""
    GaussianBeam(center, dir, waist; wavelength = 1064e-9)

A focused Gaussian beam through `center` along `dir` (normalised on construction).
"""
struct GaussianBeam
    center::Vector{Float64}
    dir::Vector{Float64}
    waist::Float64
    zR::Float64
end

GaussianBeam(center, dir, waist; wavelength = 1064e-9) =
    GaussianBeam(collect(Float64, center), normalize(collect(Float64, dir)),
                 waist, pi * waist^2 / wavelength)

"Intensity of `b` at (x, y, z) relative to its peak intensity."
@inline function beam_profile(b::GaussianBeam, x, y, z)
    dx, dy, dz = x - b.center[1], y - b.center[2], z - b.center[3]
    zax = dx * b.dir[1] + dy * b.dir[2] + dz * b.dir[3]
    r2  = (dx - zax * b.dir[1])^2 + (dy - zax * b.dir[2])^2 + (dz - zax * b.dir[3])^2
    w2  = b.waist^2 * (1 + (zax / b.zR)^2)
    return (b.waist^2 / w2) * exp(-2 * r2 / w2)
end

"""
    trap_intensity(beams, peaks)

Return `(x, y, z) -> I` for the incoherent sum of `beams` with peak intensities `peaks`.
"""
function trap_intensity(beams::AbstractVector{GaussianBeam}, peaks::AbstractVector{<:Real})
    length(beams) == length(peaks) || error("one peak intensity per beam")
    bs, ps = collect(beams), collect(Float64, peaks)
    return (x, y, z) -> sum(ps[i] * beam_profile(bs[i], x, y, z) for i in eachindex(bs))
end

"""
    trap_potential(beams, peaks; α = KRb_polarizability_1064)

Return `(x, y, z) -> U` (J), the dipole potential -α I of the beams.
"""
function trap_potential(beams, peaks; α = KRb_polarizability_1064)
    I = trap_intensity(beams, peaks)
    return (x, y, z) -> -α * I(x, y, z)
end

# ---------------------------------------------------------------------------
# Thermal sampling
# ---------------------------------------------------------------------------

"""
    sample_thermal(potential, mass, T, n; box, rng = Random.default_rng())

Draw `n` particles from a thermal distribution at temperature `T` in `potential`.
Positions are rejection-sampled uniformly within ±`box` (a 3-vector of half-widths,
m) with Boltzmann acceptance; velocities are Maxwell-Boltzmann.

Returns `(positions, velocities)`, each a 3×n matrix.
"""
function sample_thermal(potential, mass, T, n; box, rng = Random.default_rng())
    U0  = potential(0.0, 0.0, 0.0)
    pos = Matrix{Float64}(undef, 3, n)
    k = 0
    while k < n
        x = (2 * rand(rng) - 1) * box[1]
        y = (2 * rand(rng) - 1) * box[2]
        z = (2 * rand(rng) - 1) * box[3]
        if exp(-(potential(x, y, z) - U0) / (kB * T)) > rand(rng)
            k += 1
            pos[1, k], pos[2, k], pos[3, k] = x, y, z
        end
    end
    vel = randn(rng, 3, n) .* sqrt(kB * T / mass)
    return pos, vel
end

"""
    thermal_stats(pos, vel, mass)

Kinetic temperature (K) along x, y, z, and the harmonic trap frequencies (Hz, not
angular) implied by the position spread at that temperature.
"""
function thermal_stats(pos, vel, mass)
    T = [mean(vel[i, :] .^ 2) * mass / kB for i in 1:3]
    f = [1 / std(pos[i, :] ./ sqrt(T[i] * kB / mass)) / (2pi) for i in 1:3]
    return (temperature = T, trap_freq = f)
end

# ---------------------------------------------------------------------------
# Microwave drive
# ---------------------------------------------------------------------------

const _pol = (pi          = ComplexF64[0, 0, 1],
              sigma_plus  = -1 / sqrt(2) * ComplexF64[1, 1im, 0],
              sigma_minus =  1 / sqrt(2) * ComplexF64[1, -1im, 0])

"Spherical components (π, σ+, σ-) of a Cartesian polarisation vector."
spherical_components(E) = ComplexF64[dot(_pol.pi, E), dot(_pol.sigma_plus, E), dot(_pol.sigma_minus, E)]

"Rotation matrix about `axis` by `theta` (rad), in the quaternion convention used historically here."
function rotation_matrix(axis, theta)
    axis = axis / norm(axis)
    a = cos(theta / 2)
    b, c, d = -axis * sin(theta / 2)
    return [a*a+b*b-c*c-d*d  2*(b*c-a*d)      2*(b*d+a*c)
            2*(b*c+a*d)      a*a+c*c-b*b-d*d  2*(c*d-a*b)
            2*(b*d-a*c)      2*(c*d+a*b)      a*a+d*d-b*b-c*c]
end

"""
    MicrowaveDrive(; angle_deg, pi_time, ref_dipole = 0.33Debye,
                     calibrate_on = :sigma_minus, axis = [1,0,0], E0 = [0,0,1])

Microwave field vector, in units where `E ⋅ d` (d in C·m) is a frequency in Hz.

The polarisation is `E0` rotated by `angle_deg` about `axis`. The amplitude is set so
that a transition of dipole moment `ref_dipole`, driven by the `calibrate_on`
spherical component (`:pi`, `:sigma_plus` or `:sigma_minus`), has a π-pulse time of
`pi_time`. The coupling entering the Hamiltonian is `E ⋅ d / 2`.
"""
struct MicrowaveDrive
    E::Vector{Float64}
end

function MicrowaveDrive(; angle_deg, pi_time, ref_dipole = 0.33 * Debye,
                        calibrate_on::Symbol = :sigma_minus,
                        axis = [1.0, 0.0, 0.0], E0 = [0.0, 0.0, 1.0])
    ehat = normalize(rotation_matrix(axis, angle_deg * pi / 180) * E0)
    comp = abs(spherical_components(ehat)[(pi = 1, sigma_plus = 2, sigma_minus = 3)[calibrate_on]])
    rabi = 1 / (2 * pi_time)
    return MicrowaveDrive(rabi / (ref_dipole * comp) * ehat)
end

# ---------------------------------------------------------------------------
# Molecule at one trap intensity
# ---------------------------------------------------------------------------

"""
    manifolds(Hmol)

Index ranges `(ground, excited)` of the N=0 and N=1 states, read from N(N+1) in
the diatomic basis. The eigenbasis is energy-sorted, so the same ranges label the
N=0 and N=1 eigenstates as long as the fields do not mix rotational levels.
"""
function manifolds(Hmol)
    N2 = real.(diag(Hmol.MolOp.N[1:3] * Hmol.MolOp.N[1:3]))
    ground  = findall(x -> isapprox(x, 0; atol = 1e-8), N2)
    excited = findall(x -> isapprox(x, 2), N2)
    (ground == 1:length(ground) && excited == length(ground) .+ (1:length(excited))) ||
        error("N=0 and N=1 states are not contiguous at the start of the basis")
    return 1:length(ground), excited[1]:excited[end]
end

"""
    DressedMolecule(sol, dipOp, drive, ground, excited)

The N<=1 part of the molecule at one trap intensity, in its eigenbasis.

- `sol`     a `diatomic_jl.solve.sol` (eigenvalues in Hz, eigenvectors)
- `dipOp`   Cartesian dipole operators, `Hamiltonian.getDipoleMatrix(Hmol.MolOp)`
- `drive`   a `MicrowaveDrive`
- `ground`, `excited`  from `manifolds(Hmol)`

`V = E ⋅ d / 2` is kept in full, so it includes the small induced-dipole terms within
each manifold as well as the N=0 <-> N=1 couplings.
"""
struct DressedMolecule
    E::Vector{Float64}
    V::Matrix{ComplexF64}
    ground::UnitRange{Int}
    excited::UnitRange{Int}
    intensity::Float64
end

function DressedMolecule(sol, dipOp, drive::MicrowaveDrive, ground, excited)
    keep = 1:last(excited)
    U = sol.vec[:, keep]
    V = sum(drive.E[k] * (U' * (dipOp[k] * U)) for k in 1:3) / 2
    return DressedMolecule(sol.val[keep], Matrix(Hermitian((V + V') / 2)), ground, excited,
                           Float64(sol.Intensity[1]))
end

"""
    rotating_frame_hamiltonian(m, f)

Hamiltonian (Hz) in the frame rotating at microwave frequency `f`: N=1 energies are
shifted down by `f`, coupling `m.V`.
"""
function rotating_frame_hamiltonian(m::DressedMolecule, f)
    H = copy(m.V)
    for j in eachindex(m.E)
        H[j, j] += m.E[j] - (j in m.excited ? f : 0.0)
    end
    return Hermitian(H)
end

"Propagator exp(-2πi H t) of a Hermitian H given in Hz."
function propagator(H::Hermitian, t)
    F = eigen(H)
    return F.vectors * Diagonal(cis.(-2pi * t .* F.values)) * F.vectors'
end

_check_states(m, states, weights) =
    (length(states) == length(weights) || error("one weight per state");
     all(in(m.ground), states) || error("initial states must be in the N=0 manifold"))

"""
    spectrum(m, freqs, states, weights, t)

Excited-state (N=1) population after a square microwave pulse of length `t`, for each
initial state in `states` (N=0 indices) and each frequency in `freqs`. Uses the
whole N<=1 space. `weights` are initial populations.

Returns a `length(states) × length(freqs)` matrix.
"""
function spectrum(m::DressedMolecule, freqs, states, weights, t)
    _check_states(m, states, weights)
    out = zeros(length(states), length(freqs))
    for (i, f) in enumerate(freqs)
        U = propagator(rotating_frame_hamiltonian(m, f), t)
        for (k, g) in enumerate(states)
            out[k, i] = weights[k] * sum(abs2, @view U[m.excited, g])
        end
    end
    return out
end

"""
    spectrum_truncated(m, freqs, states, weights, t; window = 20e3)

As `spectrum`, but for each initial state `g` and frequency only `g` and the N=1
states within `window` (Hz) of resonance with it are kept. Much faster, and exact
whenever the omitted states are far detuned.
"""
function spectrum_truncated(m::DressedMolecule, freqs, states, weights, t; window = 20e3)
    _check_states(m, states, weights)
    out = zeros(length(states), length(freqs))
    Vd  = real.(diag(m.V))
    sub = Int[]
    for (i, f) in enumerate(freqs), (k, g) in enumerate(states)
        empty!(sub)
        push!(sub, g)
        for j in m.excited
            abs(m.E[j] - m.E[g] - f + Vd[j]) < window && push!(sub, j)
        end
        length(sub) == 1 && continue
        H = m.V[sub, sub]
        for a in eachindex(sub)
            H[a, a] += m.E[sub[a]] - m.E[g] - (a == 1 ? 0.0 : f)
        end
        U = propagator(Hermitian(H), t)
        out[k, i] = weights[k] * sum(abs2, @view U[2:end, 1])
    end
    return out
end

"""
    intensity_averaged_spectrum(intensities, scan, dipOp, drive, ground, excited,
                                freqs, states, weights, t; window = 20e3)

Truncated spectrum summed over molecules, each at its own trap intensity.

Each entry of `intensities` (W/m^2) is assigned to the nearest point of `scan`
(a vector of `solve.sol` from `solve.scanIntensity`); the spectrum is computed once
per occupied grid point and weighted by how many molecules fall there.

Returns `(total, counts)`: `total` is the spectrum summed over molecules and initial
states (length(freqs)); `counts[j]` is the number of molecules on `scan[j]`.
"""
function intensity_averaged_spectrum(intensities, scan, dipOp, drive::MicrowaveDrive,
                                     ground, excited, freqs, states, weights, t;
                                     window = 20e3)
    grid   = [s.Intensity[1] for s in scan]
    counts = zeros(Int, length(grid))
    for I in intensities
        counts[argmin(abs.(grid .- I))] += 1
    end
    used  = findall(>(0), counts)
    parts = zeros(length(freqs), length(used))
    @threads for u in eachindex(used)
        m = DressedMolecule(scan[used[u]], dipOp, drive, ground, excited)
        parts[:, u] = vec(sum(spectrum_truncated(m, freqs, states, weights, t; window), dims = 1))
    end
    return parts * counts[used], counts
end

# ---------------------------------------------------------------------------
# Repeated pulse + partial measurement
# ---------------------------------------------------------------------------

"""
    multistep_spectrum(m, freqs, states, weights, step_time, n_steps, fidelity;
                       first_pulse = nothing)

Accumulated detection signal for `n_steps` cycles of (microwave pulse of `step_time`,
then a measurement that detects a fraction `fidelity` of the N=1 population and
removes it from the state). If `first_pulse` is given, one extra cycle with a pulse
of that length is run first.

Returns a `length(states) × length(freqs)` matrix.
"""
function multistep_spectrum(m::DressedMolecule, freqs, states, weights,
                            step_time, n_steps, fidelity; first_pulse = nothing)
    _check_states(m, states, weights)
    out    = zeros(length(states), length(freqs))
    keep   = sqrt(1 - fidelity)
    pulses = first_pulse === nothing ? fill(step_time, n_steps) : [first_pulse; fill(step_time, n_steps)]
    for (i, f) in enumerate(freqs)
        H  = rotating_frame_hamiltonian(m, f)
        Us = Dict(t => propagator(H, t) for t in unique(pulses))
        for (k, g) in enumerate(states)
            psi = zeros(ComplexF64, length(m.E))
            psi[g] = 1
            signal = 0.0
            for t in pulses
                psi = Us[t] * psi
                signal += fidelity * sum(abs2, @view psi[m.excited])
                psi[m.excited] .*= keep
            end
            out[k, i] = weights[k] * signal
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Detection of reaction products by REMPI outside the dark line
# ---------------------------------------------------------------------------

"""
    product_speeds(N_KRb, F_Rb; exergonic = 0.05)

Speeds (m/s) of the KRb and Rb fragments when the reaction releases `exergonic` GHz,
less the rotational energy of KRb in `N_KRb` and 6.628 GHz if Rb is left in F=2.
"""
function product_speeds(N_KRb, F_Rb; exergonic = 0.05)
    B_KRb = 2.2227 / 2
    U = (F_Rb == 2 ? 6.628 : 0.0) + B_KRb * N_KRb * (N_KRb + 1)
    E = (exergonic - U) * h * 1e9
    mKRb, mRb = (40 + 87) * amu, 87 * amu
    return (KRb = sqrt(2 * mRb / (mKRb * (mKRb + mRb)) * E),
            Rb  = sqrt(2 * mKRb / (mRb * (mKRb + mRb)) * E))
end

"Geometry of the ionisation beam: waist, dark-line half-width (m), peak ionisation probability."
struct REMPIConfig
    beam_waist::Float64
    dark_line_radius::Float64
    max_ionization_prob::Float64
end

"""
    linear_detection(cfg)

Detection model `r -> (ionised fraction, surviving fraction)` where the ionisation
probability follows the beam intensity outside the dark line, and everything that is
not ionised survives.
"""
function linear_detection(cfg::REMPIConfig)
    return function (r)
        eta = r >= cfg.dark_line_radius ?
              cfg.max_ionization_prob * exp(-2 * r^2 / cfg.beam_waist^2) : 0.0
        return eta, 1 - eta
    end
end

_draw(v::Real, rng) = v
_draw(f, rng) = f(rng)

"""
    product_spectrum(m, freqs, states, weights, detection;
                     t_dark, microwave_time, step_time, n_steps, speed,
                     n_particles = 100, rng = Random.default_rng())

Monte Carlo signal from reaction products that fly out of the dark line while being
probed. Each particle is born at a uniformly random time in the last `t_dark` before
probing starts, moves radially at `speed` (a number, or a function `rng -> speed`)
with an isotropic direction projected onto the detection plane, and is probed every
`step_time` for `n_steps` cycles. Each cycle is a microwave pulse of `microwave_time`
followed by detection, where `detection(r)` returns the ionised and surviving
fractions of the N=1 population at radius `r` (e.g. `linear_detection(cfg)`).

Returns the mean signal per particle, `length(states) × length(freqs)`.
"""
function product_spectrum(m::DressedMolecule, freqs, states, weights, detection;
                          t_dark, microwave_time, step_time, n_steps, speed,
                          n_particles = 100, rng = Random.default_rng())
    _check_states(m, states, weights)

    ion  = zeros(n_particles, n_steps)
    surv = ones(n_particles, n_steps)
    for p in 1:n_particles
        v = _draw(speed, rng)
        t_spawn = rand(rng) * t_dark
        r_scale = v * sin(acos(2 * rand(rng) - 1))
        for s in 1:n_steps
            ion[p, s], surv[p, s] = detection(r_scale * ((t_dark - t_spawn) + (s - 1) * step_time))
        end
    end
    keep = sqrt.(surv)

    n   = length(m.E)
    out = zeros(length(states), length(freqs))
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @threads for i in eachindex(freqs)
            U    = propagator(rotating_frame_hamiltonian(m, freqs[i]), microwave_time)
            psi  = zeros(ComplexF64, n, n_particles)
            next = similar(psi)
            for (k, g) in enumerate(states)
                fill!(psi, 0)
                psi[g, :] .= 1
                signal = 0.0
                for s in 1:n_steps
                    mul!(next, U, psi)
                    for p in 1:n_particles
                        signal += ion[p, s] * sum(abs2, @view next[m.excited, p])
                        @views next[m.excited, p] .*= keep[p, s]
                    end
                    psi, next = next, psi
                end
                out[k, i] = weights[k] * signal / n_particles
            end
        end
    finally
        BLAS.set_num_threads(old_blas)
    end
    return out
end

end # module
