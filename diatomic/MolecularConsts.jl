
using PhysicalConstants.CODATA2018
# Constants
h = PlanckConstant.val
muN = 5.0507837461e-27
bohr = BohrRadius.val
eps0 = VacuumElectricPermittivity.val
c_l = SpeedOfLightInVacuum.val

DebyeSI = 3.33564e-30

# Rb87Cs133 Constants

Rb87Cs133 = Dict(
    "Name" => "Rb87Cs133",
    "I1" => 1.5,
    "I2" => 3.5,
    "d0" => 1.225 * DebyeSI,
    "binding" => 114268135.25e6 * h,
    "Brot" => 490.173994326310e6 * h,
    "Drot" => 207.3 * h,
    "Q1" => -809.29e3 * h,
    "Q2" => 59.98e3 * h,
    "C1" => 98.4 * h,
    "C2" => 194.2 * h,
    "C3" => 192.4 * h,
    "C4" => 19.0189557e3 * h,
    "MuN" => 0.0062 * muN,
    "Mu1" => 1.8295 * muN,
    "Mu2" => 0.7331 * muN,
    "a0" => 2020*4*pi*eps0*bohr^3, #1064nm
    "a2" => 1997*4*pi*eps0*bohr^3, #1064nm
    "Beta" => 0
)
# K41Cs133 Constants
K41Cs133 = Dict(
    "Name" => "K41Cs133",
    "I1" => 1.5,
    "I2" => 3.5,
    "d0" => 1.84 * DebyeSI,
    "Brot" => 880.326e6 * h,
    "Drot" => 0 * h,
    "Q1" => -0.221e6 * h,
    "Q2" => 0.075e6 * h,
    "C1" => 4.5 * h,
    "C2" => 370.8 * h,
    "C3" => 9.9 * h,
    "C4" => 628 * h,
    "MuN" => 0.0 * muN,
    "Mu1" => 0.143 * (1 - 1340.7e-6) * muN,
    "Mu2" => 0.738 * (1 - 6337.1e-6) * muN,
    "a0" => 7.783e6 * h,
    "a2" => 0,
    "Beta" => 0
)

# K40Rb87 Constants
K40Rb87 = Dict(
    "Name" => "K40Rb87",
    "I1" => 4,
    "I2" => 1.5,
    "d0" => 0.566 * DebyeSI,
    "Brot" => 1113.950e6 * h,
    "Drot" => 0 * h,
    "Q1" => 0.45e6 * h,
    "Q2" => -1.41e6 * h,
    "C1" => -24.1 * h,
    "C2" => 419.5 * h,
    "C3" => -48.2 * h,
    "C4" => -2028.8 * h,
    "MuN" => 0.0140 * muN,
    "Mu1" => -0.324 * (1 - 1321e-6) * muN,
    "Mu2" => 1.834 * (1 - 3469e-6) * muN,
    "a0" => 55.3*1e-4*h*eps0*c_l*2, #0,  #1/3*(100 +2*33) *h *1e-4 , #5.53e-5 * 1e6 * h,
    "a2" => 44.7*1e-4*h*eps0*c_l*2, #2/3*(100 - 33)* h *1e-4*eps0*c_l/2, #4.47e-5 * 1e6 * h,
    "Beta" => 0
)


K40Rb87T = Dict(
    "Name" => "K40Rb87_Till",
    "I1" => 4,
    "I2" => 1.5,
    "d0" => 0.573999 * DebyeSI,
    "Brot" => 1113.9514e6 * h,
    "Drot" => 0 * h,
    "Q1" => 0.452e6 * h,
    "Q2" => -1.308e6 * h,
    "C1" => -24.1 * h,
    "C2" => 420.1 * h,
    "C3" => -48.2 * h,
    "C4" => -2030.4 * h,
    "MuN" => 0.0140 * muN,
    "Mu1" => -0.324 * (1 - 1321e-6) * muN,
    "Mu2" => 1.834 * (1 - 3469e-6) * muN,
    "a0" => 0, #55.3*1e-4*h*eps0*c_l*2, #0,  #1/3*(100 +2*33) *h *1e-4 , #5.53e-5 * 1e6 * h,
    "a2" => 44.7*1e-4*h*eps0*c_l*2, #2/3*(100 - 33)* h *1e-4*eps0*c_l/2, #4.47e-5 * 1e6 * h,
    "Beta" => 0
)


# Na23K40 Constants
Na23K40 = Dict(
    "Name" => "Na23K40",
    "I1" => 1.5,
    "I2" => 4,
    "d0" => 2.72 * DebyeSI,
    "Brot" => 2.8217297e9 * h,
    "Drot" => 0 * h,
    "Q1" => -0.187e6 * h,
    "Q2" => 0.899e6 * h,
    "C1" => 117.4 * h,
    "C2" => -97.0 * h,
    "C3" => -48.4 * h,
    "C4" => -409 * h,
    "MuN" => 0.0253 * muN,
    "Mu1" => 1.477 * (1 - 624.4e-6) * muN,
    "Mu2" => -0.324 * (1 - 1297.4e-6) * muN,
    "a0" => 0 * h,
    "a2" => 0 * h,
    "Beta" => 0
)

# Na23Rb87 Constants
Na23Rb87 = Dict(
    "Name" => "Na23Rb87",
    "I1" => 1.5,
    "I2" => 1.5,
    "d0" => 3.2 * DebyeSI,
    "Brot" => 2.0896628e9 * h,
    "Drot" => 0 * h,
    "Q1" => -0.139e6 * h,
    "Q2" => -3.048e6 * h,
    "C1" => 60.7 * h,
    "C2" => 983.8 * h,
    "C3" => 259.3 * h,
    "C4" => 6.56e3 * h,
    "MuN" => 0.001 * muN,
    "Mu1" => 1.484 * muN,
    "Mu2" => 1.832 * muN,
    "a0" => 0 * h,
    "a2" => 0 * h,
    "Beta" => 0
)


