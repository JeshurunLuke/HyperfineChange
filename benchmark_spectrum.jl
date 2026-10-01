# Benchmark script comparing original vs optimized simulateSpectrumWrap
using BenchmarkTools

# Include the optimized version
include("IntDist_optimized.jl")

# Initialize constants (you'll need to have these functions available)
# SpectrumConstants.init_constants!(rotation_matrix, axisX, Ecomp1, decompose_spherical)

"""
Run this after loading your notebook to compare performance:

# Setup (from your notebook)
mFieldScan = [2228.08...]*1e6
statePop = [8]
popValues = ones(size(statePop))
driveTime = 200e-6

# Get the necessary data for a single intensity point
keyInd = 1
dipOpU = UnitaryDict[collect(keys(UnitaryDict))[keyInd]].dipOpU
H0 = Diagonal(IntensityScan[keyInd].val)
N_product = Hmol.MolOp.N[1:3] * Hmol.MolOp.N[1:3]
Hmicr = zeros(Float64, size(N_product))
Hmicr[findall(x->isapprox(x, 2), N_product)] .= 1

# Initialize constants
include("IntDist_optimized.jl")
SpectrumConstants.init_constants!(rotation_matrix, axisX, Ecomp1, decompose_spherical)

# Benchmark original
println("Original version:")
@btime simulateSpectrumWrap(\$mFieldScan, \$statePop, \$popValues, \$driveTime,
                            \$Hmicr, \$H0, \$dipOpU; truncateEnergy = 20e3)

# Benchmark optimized
println("\\nOptimized version:")
@btime simulateSpectrumWrap_optimized(\$mFieldScan, \$statePop, \$popValues, \$driveTime,
                                      \$Hmicr, \$H0, \$dipOpU; truncateEnergy = 20e3)

# Check results are the same
result_original = simulateSpectrumWrap(mFieldScan, statePop, popValues, driveTime,
                                       Hmicr, H0, dipOpU; truncateEnergy = 20e3)
result_optimized = simulateSpectrumWrap_optimized(mFieldScan, statePop, popValues, driveTime,
                                                   Hmicr, H0, dipOpU; truncateEnergy = 20e3)

println("\\nResults match: ", isapprox(result_original, result_optimized, rtol=1e-10))
println("Max difference: ", maximum(abs.(result_original .- result_optimized)))

# Profile optimized version
using Profile, PProf
@profile simulateSpectrumWrap_optimized(mFieldScan, statePop, popValues, driveTime,
                                        Hmicr, H0, dipOpU; truncateEnergy = 20e3)
pprof()
"""

# Allocation comparison function
function compare_allocations(mFieldScan, statePop, popValues, driveTime, Hmicr, H0, dipOpU)
    println("=== Allocation Comparison ===\n")

    println("Original version:")
    stats_orig = @timed simulateSpectrumWrap(mFieldScan, statePop, popValues, driveTime,
                                             Hmicr, H0, dipOpU; truncateEnergy = 20e3)
    println("  Time: $(stats_orig.time) s")
    println("  Allocations: $(stats_orig.bytes / 1e6) MB")
    println("  GC time: $(stats_orig.gctime) s")

    println("\nOptimized version:")
    stats_opt = @timed simulateSpectrumWrap_optimized(mFieldScan, statePop, popValues, driveTime,
                                                      Hmicr, H0, dipOpU; truncateEnergy = 20e3)
    println("  Time: $(stats_opt.time) s")
    println("  Allocations: $(stats_opt.bytes / 1e6) MB")
    println("  GC time: $(stats_opt.gctime) s")

    println("\nSpeedup: $(stats_orig.time / stats_opt.time)x")
    println("Memory reduction: $(stats_orig.bytes / stats_opt.bytes)x")
end
