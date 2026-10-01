# Optimized version of simulateSpectrumWrap
using LinearAlgebra

# Precompute constants outside the function
module SpectrumConstants
    using LinearAlgebra

    const DebyeSI = 3.33564e-30
    const TRUNCATE = 144
    const ROTATION_ANGLE = 87.5
    const PI_TIME = 398e-6
    const RABI_RATE = 1/(2*PI_TIME)

    # These should be set once based on your beam geometry
    # You'll need to call init_constants! before using the optimized function
    E_field_res::Vector{Float64} = Float64[]
    twopi::Float64 = 2π
    inv_twopi::Float64 = 1/(2π)

    function init_constants!(rotation_matrix_func, axisX, Ecomp1, decompose_spherical_func)
        Erot = rotation_matrix_func(axisX, ROTATION_ANGLE*π/180) * Ecomp1
        Ecomp = normalize!(Erot)
        sigmaPcomp = abs(decompose_spherical_func(Erot)[3])
        E_field_hz = RABI_RATE/(0.33*DebyeSI)*(1/sigmaPcomp)

        resize!(E_field_res, 3)
        E_field_res .= E_field_hz .* Ecomp
        nothing
    end
end


function simulateSpectrumWrap_optimized(
    mField::AbstractVector,
    statePop::AbstractVector{Int},
    popValues::AbstractVector,
    driveTime::Float64,
    Hmicr::Matrix{Float64},
    H0::Diagonal{Float64},
    dipOpU::Vector;
    truncateEnergy::Float64 = 20e3
)
    # Pre-allocate all arrays we'll need
    n_mField = length(mField)
    n_states = size(H0, 1)

    popN1 = zeros(Float64, 36, n_mField)

    # Preallocate working arrays
    H_full = zeros(ComplexF64, n_states, n_states)
    H_diag = zeros(Float64, n_states)
    dipole_sum = zeros(ComplexF64, n_states, n_states)

    # Precompute dipole operator sum (this is constant across all iterations!)
    E_field_res = SpectrumConstants.E_field_res
    @inbounds for k in 1:3
        @. dipole_sum += E_field_res[k] * dipOpU[k]
    end

    # Constant for truncation
    const statesOI_range = 1:144
    const excited_threshold = 36

    # Main computation loop
    @inbounds for (i, m_i) in enumerate(mField)
        for (_ind, stateGS) in enumerate(statePop)
            # Build Hamiltonian in-place
            # H = 2π*(H0 - m_i*Hmicr + dipole_sum - Diagonal(H0[stateGS, stateGS]))

            H0_shift = H0.diag[stateGS]

            # Construct H efficiently
            @inbounds for j in 1:n_states
                for k in 1:n_states
                    H_full[j, k] = SpectrumConstants.twopi * (
                        (j == k ? H0.diag[j] - H0_shift : 0.0)
                        - m_i * Hmicr[j, k]
                        + dipole_sum[j, k]
                    )
                end
                H_diag[j] = real(H_full[j, j]) * SpectrumConstants.inv_twopi
            end

            # Find states within energy threshold
            # More efficient than findall + intersect
            statesTruncated = Int[]
            sizehint!(statesTruncated, 100)  # Pre-allocate capacity

            # Always include ground state
            push!(statesTruncated, stateGS)

            # Add excited states within energy window
            for state in statesOI_range
                if state > excited_threshold && abs(H_diag[state]) < truncateEnergy
                    push!(statesTruncated, state)
                end
            end

            n_trunc = length(statesTruncated)

            # Extract truncated Hamiltonian (unavoidable allocation, but minimized)
            Htrunc = zeros(ComplexF64, n_trunc, n_trunc)
            @inbounds for (j_idx, j) in enumerate(statesTruncated)
                for (k_idx, k) in enumerate(statesTruncated)
                    Htrunc[j_idx, k_idx] = H_full[j, k]
                end
            end

            # Matrix exponential (expensive but necessary)
            Uprop = exp(-1im * Htrunc * driveTime)

            # Initial state (ground state is first in truncated basis)
            stateInit = zeros(ComplexF64, n_trunc)
            stateInit[1] = popValues[_ind]

            # Propagate
            stateFinal = Uprop' * stateInit

            # Calculate population in excited states (skip first = ground state)
            pop_excited = 0.0
            @inbounds for j in 2:n_trunc
                pop_excited += abs2(stateFinal[j])
            end

            popN1[stateGS, i] = pop_excited
        end
    end

    # Return sum over initial states
    return sum(popN1, dims=1)
end


# Version with even more aggressive optimizations for single magnetic field point
function simulateSpectrumWrap_single_field(
    m_field::Float64,
    stateGS::Int,
    popValue::Float64,
    driveTime::Float64,
    Hmicr::Matrix{Float64},
    H0::Diagonal{Float64},
    dipOpU::Vector;
    truncateEnergy::Float64 = 20e3
)
    n_states = size(H0, 1)

    # Preallocate
    H_full = zeros(ComplexF64, n_states, n_states)
    dipole_sum = zeros(ComplexF64, n_states, n_states)

    # Precompute dipole sum
    E_field_res = SpectrumConstants.E_field_res
    @inbounds for k in 1:3
        @. dipole_sum += E_field_res[k] * dipOpU[k]
    end

    # Build Hamiltonian
    H0_shift = H0.diag[stateGS]

    @inbounds for j in 1:n_states
        for k in 1:n_states
            H_full[j, k] = SpectrumConstants.twopi * (
                (j == k ? H0.diag[j] - H0_shift : 0.0)
                - m_field * Hmicr[j, k]
                + dipole_sum[j, k]
            )
        end
    end

    # Truncate to relevant states
    statesTruncated = Int[stateGS]
    @inbounds for state in 37:144  # Excited states only
        if abs(real(H_full[state, state]) * SpectrumConstants.inv_twopi) < truncateEnergy
            push!(statesTruncated, state)
        end
    end

    n_trunc = length(statesTruncated)
    Htrunc = zeros(ComplexF64, n_trunc, n_trunc)
    @inbounds for (j_idx, j) in enumerate(statesTruncated)
        for (k_idx, k) in enumerate(statesTruncated)
            Htrunc[j_idx, k_idx] = H_full[j, k]
        end
    end

    Uprop = exp(-1im * Htrunc * driveTime)

    stateInit = zeros(ComplexF64, n_trunc)
    stateInit[1] = popValue

    stateFinal = Uprop' * stateInit

    # Sum excited state population
    pop_excited = 0.0
    @inbounds for j in 2:n_trunc
        pop_excited += abs2(stateFinal[j])
    end

    return pop_excited
end


# Wrapper for integration with existing code
function wrapIntensity_optimized(
    IntensityDist,
    IntensityScan,
    UnitaryDict,
    Hmol,
    mField,
    statePop,
    popValues,
    driveTime;
    truncateEnergy = 20e3
)
    keysOI = collect(keys(UnitaryDict))
    res = zeros(Float64, length(IntensityDist), length(mField))

    N_product = Hmol.MolOp.N[1:3] * Hmol.MolOp.N[1:3]
    Hmicr = zeros(Float64, size(N_product))
    Hmicr[findall(x->isapprox(x, 2), N_product)] .= 1

    # Pre-compute key indices to avoid repeated searches
    key_indices = [argmin(abs.(keysOI .- intensity))[1] for intensity in IntensityDist]

    for (indC, keyInd) in enumerate(key_indices)
        dipOpU = UnitaryDict[keysOI[keyInd]].dipOpU
        H0 = Diagonal(IntensityScan[keyInd].val)

        res[indC, :] .= simulateSpectrumWrap_optimized(
            mField, statePop, popValues, driveTime,
            Hmicr, H0, dipOpU;
            truncateEnergy = truncateEnergy
        )
    end

    return res
end
