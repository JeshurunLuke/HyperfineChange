module diatomic_jl


include("QuantumToolbox/QuantumWrapper.jl")
using .QuantumWrapper
import .QuantumWrapper: AM_Toolbox
import .QuantumWrapper.State_Toolbox.State
import .QuantumWrapper.State_Toolbox.KetName

#using Reexport
#@reexport using .QuantumWrappera
export Molecule, Hamiltonian, solve, plotting, calculate, AM_Toolbox, getBasisUC, State, KetName

function getBasisUC end



module MoleculeTypes
include("diatomic/MolecularConsts.jl")
export Rb87Cs133, K41Cs133, K40Rb87 ,Na23K40, Na23Rb87, K40Rb87T

end



### Initializes the molecule and the Appropriate spin operators
module Molecule
import ..QuantumWrapper.AM_Toolbox: Basis, Node, angularMomentum, couple, endNode, AngularNode
using ..QuantumWrapper.SpinOperator
using LinearAlgebra
using SparseArrays
export moleculeProperties, generatorMolecule
struct moleculeProperties{T<:Node}
    Constants
    basisTree::T
    N
    I1
    I2
end

function generateMolecule(MoleculeInfo, Nmax)
    I1 = angularMomentum("I1", [MoleculeInfo["I1"]])
    I2 = angularMomentum("I2", [ MoleculeInfo["I2"]])
    N = angularMomentum("N", [0:Nmax...])
    I_c = couple("I_c", I2, I1)
    F = couple("F", N, I_c)
    Nop, I2op, I1op = generateSpinOperator(F)
    moleculeProperties(MoleculeInfo,F, Nop, I1op, I2op)
end
function generateSpinOperator(F::Node)
    EndNodes = endNode(F)
    Operators = []
    for (i, node) in enumerate(EndNodes)
        dimsA = i == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(i+ 1):end]])
        dimsB = i == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(i-1):-1:1]])



        I_a = sparse(1.0I, dimsA, dimsA)
        I_b = sparse(1.0I, dimsB, dimsB)

        push!(Operators, [kron(I_b, kron(x_operator(node), I_a)), kron(I_b, kron(y_operator(node), I_a)), kron(I_b, kron(z_operator(node), I_a)), kron(I_b, kron(raising_operator(node), I_a))])
    end
    Operators
end
end
##### Used to Generate the Hamiltonian with the spin operators ################################
module Hamiltonian
import ..QuantumWrapper.AM_Toolbox: endNode, getBasis#, getBasisUC
import ..QuantumWrapper: AM_Toolbox
export generateHamiltonian, OpticalBeam, MoleculeHamiltonian, zeeman_ham, dc, DipoleMatrix, getDipoleMatrix
import ..Molecule: moleculeProperties
include("diatomic/Hamiltonian.jl")
using ..diatomic_jl
import ..diatomic_jl: getBasisUC
getBasisUC(mol::MoleculeHamiltonian) = AM_Toolbox.getBasisUC(mol.MolOp.basisTree)


function AC_Stark(M::moleculeProperties, Beam::OpticalBeam)
    mag = sqrt(2*Beam.Intensity/(PhysConstants["c"]*PhysConstants["e0"]))#/10000#/sqrt(2)
    a0, a2 = M.Constants["a0"], M.Constants["a2"]
    dims = prod([length(node.spin) for node in endNode(M.basisTree)])
    pol = normalize(Beam.pol_space)
    
    A2 = sqrt(3/2)*T2_C(M)
    P2 = [scaler*sparse(1.0I, dims, dims) for scaler in makeT2(pol, conj.(pol))]
    -mag^2*(sparse(1.0I, dims, dims)*a0 + tensor_dot(A2, P2, order = 2)*a2)/4
end




end

##### Used to solve the Hamiltonian ################################
module solve
export scanZeeman, scanStark, diagonalize, sol, scanIntensity
import ..Hamiltonian: MoleculeHamiltonian
using LinearAlgebra
using Arpack
using SparseArrays


mutable struct sol
    B_field
    E_field
    Intensity
    vec
    val
end

function diagonalize(H::MoleculeHamiltonian{T}, B_field, E_field; Intensity = 1,  nev = 36) where {T<:AbstractMatrix}
    J2Hz = 1.509190311676e33

    vals, vecs = eigen(Hermitian(Matrix(H.Hhfs .+ H.Hzeem*B_field .+ E_field*H.Hdc + sum(Intensity .* H.Hac))))
    sorted_indices = sortperm(vals)
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(B_field, E_field, [Intensity...], vecs, vals*J2Hz)
end
function diagonalize(H::MoleculeHamiltonian{T}, B_field, E_field; Intensity = 1,  nev = 36) where T<:SparseMatrixCSC
    J2Hz = 1.509190311676e33
    vals, vecs = eigs(H.Hhfs .+ H.Hzeem*B_field .+ E_field*H.Hdc + sum(Intensity .* H.Hac), nev = nev, which=:SM)
    sorted_indices = sortperm(real.(vals))
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(B_field, E_field, [Intensity...], vecs, real.(vals)*J2Hz)
end

function scanIntensity(H::MoleculeHamiltonian, B_Field::T, E_Field::T, Intensity) where T<:Number
    J2Hz = 1.509190311676e33
    solutions = Vector{sol}(undef, size(Intensity, 1))
    for (i, Int_i) in enumerate(Intensity)
        solutions[i] = diagonalize(H, B_Field, E_Field, Intensity = Int_i)
    end
    solutions
end


function scanZeeman(H::MoleculeHamiltonian, B_field::Vector{<:Real}, E_field::T; Intensity = 1, nev = 36) where T<:Number
    J2Hz = 1.509190311676e33
    solutions = Vector{sol}(undef, length(B_field))
    for (i, B_i) in enumerate(B_field)
        solutions[i] = diagonalize(H, B_i , E_field, Intensity = 1, nev = nev)
    end
    solutions
end
function scanStark(H::MoleculeHamiltonian, B_field::T, E_field::Vector{<:Real}; Intensity = 1, nev = 36) where T<:Number
    J2Hz = 1.509190311676e33
    solutions = Vector{sol}(undef, length(E_field))
    for (i, E_i) in enumerate(E_field)
        solutions[i] = diagonalize(H, B_field , E_i, Intensity= 1, nev = nev)
    end
    solutions
end
end

######## Plotting Utils ###############


module calculate
export findState, findTransition, transition_dipole_moment, electric_moment, magnetic_moment, diabaticRamp, findMaxOverlap, findAvoidedCrossing, scanDiabiaticEnergy, diabaticOrder
using ..QuantumWrapper.State_Toolbox
import ..QuantumWrapper.AM_Toolbox: Basis, Node, getBasisUC, endNode
import ..Hamiltonian: MoleculeHamiltonian, DipoleMatrix, zeeman_ham, dc

import ..solve: sol
import ..QuantumWrapper.State_Toolbox.State
import ..QuantumWrapper.State_Toolbox.Ket
import ..QuantumWrapper.State_Toolbox.KetName
include("diatomic/calculate.jl")
findMaxOverlap(basisState::State, eigvec) = findMaxOverlap(Ket(basisState), eigvec)
findMaxOverlap(basisState::State, eigvec, ind) = findMaxOverlap(Ket(basisState), eigvec)
function findMaxOverlap(basisState::Vector{<:ComplexF64}, eigvec, indS)
    argmax(abs.(transpose(eigvec)*basisState))
end
function findMaxOverlap(basisState::Vector{<:ComplexF64}, eigvec)
    #=
    overlap = 0
    
    indOI = 0 

    for ind in 1:size(eigvec, 2)
        overlapcurr = abs2.(dot(basisState, eigvec[:, ind]))
        if overlapcurr >= overlap
            overlap = overlapcurr
            indOI = ind
        end
    end
    indOI
    =#
    argmax(abs.(transpose(eigvec)*basisState))

end


function magnetic_moment(H::MoleculeHamiltonian, eigsols::Vector{sol}; rev = true) # Returns the magnetic moment with each row as an eigenstate
    uN = 5.0507837461e-27
    sol_vec = diabaticOrder(H, eigsols, rev = rev)
    mag_mom = []
    for sol_i in sol_vec
        push!(mag_mom, real.(magnetic_moment(H, sol_i.vec))./uN)
    end
    mag_mom = hcat(mag_mom...)
end

function electric_moment(H::MoleculeHamiltonian, eigsols::Vector{sol}; rev = true)
    debye = 3.33564e-30
    sol_vec = diabaticOrder(H, eigsols, rev = rev)
    dip_mom = []
    for sol_i in sol_vec
        push!(dip_mom, real.(electric_moment(H, sol_i.vec))./debye)
    end
    dip_mom = hcat(dip_mom...) 
end



function diabaticOrder(Hmol::MoleculeHamiltonian, eigsol_vec::Vector{sol}; rev = true)
    iters = length([i for i in eigsol_vec])
    diabaticVec = sol[]
    if rev
        push!(diabaticVec, eigsol_vec[end])
        for (indOI, sol_i) in enumerate(reverse(eigsol_vec[1:(iters-1)]))
            indOrder = []
            for (start, col) in enumerate(eachcol(diabaticVec[end].vec))
                push!(indOrder, findMaxOverlap(Vector(col), sol_i.vec, start))
            end
            push!(diabaticVec,sol(sol_i.B_field, sol_i.E_field, sol_i.Intensity, sol_i.vec[:, indOrder], sol_i.val[indOrder]))
        end
        return reverse(diabaticVec)
    else 
        push!(diabaticVec, eigsol_vec[1])
        for sol_i in eigsol_vec[2:end]
            indOrder = []
            for (start, col) in enumerate(eachcol(diabaticVec[end].vec))
                push!(indOrder, findMaxOverlap(Vector(col), sol_i.vec, start))
            end
            push!(diabaticVec,sol(sol_i.B_field, sol_i.E_field, sol_i.Intensity, sol_i.vec[:, indOrder], sol_i.val[indOrder]))
        end
        return diabaticVec
    end

end


State(comp::Array{ComplexF64, 2}, QN::Array{Float64, 2}, mol::MoleculeHamiltonian) = State(normalize!(comp), QN, getBasisUC(mol.MolOp.basisTree))
State(QN::Vector{<:Float64}, mol::MoleculeHamiltonian) = State([1.0], [QN], getBasisUC(mol.MolOp.basisTree) )
State(comp::Vector{<:ComplexF64}, mol::MoleculeHamiltonian) = State(comp, getBasisUC(mol.MolOp.basisTree) , getBasisUC(mol.MolOp.basisTree) )
KetName(state::Vector{<:ComplexF64}, mol::MoleculeHamiltonian; QMorder = [5, 3, 2, 1]) = KetName(state, getBasisUC(mol.MolOp.basisTree), QMorder = QMorder )
KetName(state::State, mol::MoleculeHamiltonian; QMorder = [5, 3, 2, 1]) = KetName(Ket(state), getBasisUC(mol.MolOp.basisTree), QMorder = QMorder )
KetName(state::State) = KetName(Ket(state), state.basis)
end
module plotting
export plotZeemanMap, plotStarkMap, plotIntensityScan, plotTransitionPlot, plotMagneticMoment, plotDiabiatic
using ..QuantumWrapper.State_Toolbox
using LinearAlgebra
using ..QuantumWrapper.AM_Toolbox
import ..Hamiltonian: MoleculeHamiltonian
import ..QuantumWrapper.State_Toolbox.Ket
import ..solve: sol
using ..calculate
import CairoMakie
#import GLMakie
include("diatomic/plotting.jl")

function plotMagneticMoment(mol::MoleculeHamiltonian, sol_vec::Vector{sol}, var::String; N = [0], Beam = nothing)
    x = Vector{Float64}(undef, length([1 for i in sol_vec]))
    if var == "I"
        x = [sol_i.Intensity[Beam] for sol_i in sol_vec]
        xlabel = "Intensity (W/m)"
    elseif  var =="B"
        x = [sol_i.B_field for sol_i in sol_vec]*1e4
        xlabel = "B-Field (G)"
    elseif var == "E"
        x = [sol_i.E_field for sol_i in sol_vec]*1e-2
        xlabel = "E-Field (V/cm)"
    end


    mag_plot = magnetic_moment(mol, sol_vec)
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        addtraces!(p, scatter(x = x, y = vec(mag_plot[state, :]), text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = "Magnetic Moment", xlabel = "B-Field (G)", ylabel = "Magnetic Moment (muN)")
    p
end


function plotDiabiatic(mol::MoleculeHamiltonian, sol_vec::Vector{sol}, var::String; N = [0], Beam = nothing, rev = true)
    x = Vector{Float64}(undef, length([1 for i in sol_vec]))
    if var == "I"
        x = [sol_i.Intensity[Beam] for sol_i in sol_vec]
        xlabel = "Intensity (W/m)"
        title = "Intensity Scan"
    elseif  var =="B"
        x = [sol_i.B_field for sol_i in sol_vec]*1e4
        xlabel = "B-Field (G)"
        title = "Zeeman Scan"
    elseif var == "E"
        x = [sol_i.E_field for sol_i in sol_vec]*1e-2
        xlabel = "E-Field (V/cm)"
        title = "Stark Scan"
    end

    sol_vec = diabaticOrder(mol, sol_vec, rev = rev)

    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        addtraces!(p, scatter(x = x, y = vec([sol_i.val[state] for sol_i in sol_vec]), text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = title, xlabel = xlabel, ylabel = "Energy (GHz)")
    p
end
function plotElectricMoment(mol::MoleculeHamiltonian, sol_vec::Vector{sol}, var::String; N = [0], Beam = nothing)
    
    x = Vector{Float64}(undef, length([1 for i in sol_vec]))
    if var == "I"
        x = [sol_i.Intensity[Beam] for sol_i in sol_vec]
        xlabel = "Intensity (W/m)"
    elseif  var =="B"
        x = [sol_i.B_field for sol_i in sol_vec]*1e4
        xlabel = "B-Field (G)"
    elseif var == "E"
        x = [sol_i.E_field for sol_i in sol_vec]*1e-2
        xlabel = "E-Field (V/cm)"
    end


    mag_plot = electric_moment(mol, sol_vec)
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        addtraces!(p, scatter(x = x, y = vec(mag_plot[state, :]), text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = "Electric Moment", xlabel = xlabel, ylabel = "Electric Moment (debye)")
    p
end

plotTransitionPlot(mol::MoleculeHamiltonian, stateOI::State, sol_vec::Vector{sol}, var::String; Adiabatic = true, rev = true,  N = [0], Beam = 1, energyRef = 0 , comp = [1.0, 1.0, 1.0]) = plotTransitionPlot(mol, Ket(stateOI), sol_vec, var; Adiabatic = Adiabatic, rev = rev,  N = N, Beam = Beam, energyRef = energyRef , comp = comp)
plotTransitionPlot(mol::MoleculeHamiltonian, stateOI::State,  stateFinal::State, sol_vec::Vector{sol}, var::String; Adiabatic = true, rev = true, Beam = 1, energyRef = 0 , comp = [1.0, 1.0, 1.0]) = plotTransitionPlot(mol, Ket(stateOI), Ket(stateFinal), sol_vec, var; Adiabatic = Adiabatic, rev = rev, Beam = Beam, energyRef = energyRef , comp = comp)
function plotTransitionPlot(mol::MoleculeHamiltonian, stateOI::Vector{<:ComplexF64}, stateFinal::Vector{<:ComplexF64}, sol_vec::Vector{sol}, var::String; Adiabatic = true, rev = true, Beam = 1, energyRef = 0 , comp = [1.0, 1.0, 1.0])
    x = Vector{Float64}(undef, length([1 for i in sol_vec]))

    plotObj = CairoMakie

    fig = plotObj.Figure()
    ax = plotObj.Axis(fig[1, 1], ylabel = "Energy (GHz)")
    if var == "I"
        x = [sol_i.Intensity[Beam] for sol_i in sol_vec]
        ax.xlabel = "Intensity (W/m)"
    elseif  var =="B"
        x = [sol_i.B_field for sol_i in sol_vec]*1e-4
        ax.xlabel = "B-Field (G)"
    elseif var == "E"
        x = [sol_i.E_field for sol_i in sol_vec]*1e-2
        ax.xlabel = "E-Field (V/cm)"
    end
    x = rev ? reverse(x) : x
    startingState = rev ? sol_vec[end].vec[:, calculate.findMaxOverlap(stateOI, sol_vec[end].vec)] : sol_vec[1].vec[:,calculate.findMaxOverlap(stateOI, sol_vec[1].vec)]
    startingInd  = rev ? calculate.findMaxOverlap(stateOI, sol_vec[end].vec) : calculate.findMaxOverlap(stateOI, sol_vec[1].vec)
    transition_dipole_moment_vec = []
    mag = 0
    for (indF, field) in enumerate(x)
        stateC = Adiabatic ? (rev ? reverse(sol_vec)[indF].vec[:, startingInd] : sol_vec[indF].vec[:, startingInd])  : diabaticRamp(startingState, sol_vec, [x[1], field], Field = var)
        #println(typeof(stateC))
        #println(calculate.transition_dipole_moment(mol, stateC, rev ? reverse(sol_vec)[indF] : sol_vec[indF], comp = comp))
        tdm_field =  sum(abs.(calculate.transition_dipole_moment(mol, stateC, rev ? reverse(sol_vec)[indF].vec : sol_vec[indF].vec, comp = normalize!(comp))), dims = 2)
        if maximum(tdm_field) > mag
            mag = maximum(tdm_field)
        end
        push!(transition_dipole_moment_vec,tdm_field)
    end


    x = rev ? reverse(x) : x
    stateInds = [calculate.findMaxOverlap(stateFinal, sol_v.vec) for sol_v in sol_vec]
    tdms = [tdm[stateInds[ind]] for (ind, tdm) in enumerate(transition_dipole_moment_vec)]
    
    color = rev ? log10.(reverse(tdms) ./ mag)  : log10.(tdms ./ mag)


    lplot = plotObj.lines!(ax, x, [sol_i.val[stateInds[ind]] - energyRef for (ind, sol_i) in enumerate(sol_vec)]*1e-9; color = color, colormap = plotObj.Reverse(:tokyo25), colorrange = (log10(1e-2), log10(1)))

    custom_formatter(values) = map(
        v -> "10" * plotObj.Makie.UnicodeFun.to_superscript(round(Int64, v)),
        values
    )
    
    cb = plotObj.Colorbar(fig[1, 2], lplot,   ticks =log10.([1e-2, 1e-1, 1e0]), tickformat=custom_formatter,  label = "Transition Dipole Moment", vertical = true)
    fig
end
function plotTransitionPlot(mol::MoleculeHamiltonian, stateOI::Vector{<:ComplexF64}, sol_vec::Vector{sol}, var::String; Adiabatic = true, rev = true,  N = [0], Beam = 1, energyRef = 0 , comp = [1.0, 1.0, 1.0])
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    x = Vector{Float64}(undef, length([1 for i in sol_vec]))

    plotObj = CairoMakie

    fig = plotObj.Figure()
    ax = plotObj.Axis(fig[1, 1], ylabel = "Energy (GHz)")
    if var == "I"
        x = [sol_i.Intensity[Beam] for sol_i in sol_vec]
        ax.xlabel = "Intensity (W/m)"
    elseif  var =="B"
        x = [sol_i.B_field for sol_i in sol_vec]*1e4
        ax.xlabel = "B-Field (G)"
    elseif var == "E"
        x = [sol_i.E_field for sol_i in sol_vec]*1e-2
        ax.xlabel = "E-Field (V/cm)"
    end
    x = rev ? reverse(x) : x
    startingState = rev ? sol_vec[end].vec[:, calculate.findMaxOverlap(stateOI, sol_vec[end].vec)] : sol_vec[1].vec[:,calculate.findMaxOverlap(stateOI, sol_vec[1].vec)]
    startingInd  = rev ? calculate.findMaxOverlap(stateOI, sol_vec[end].vec) : calculate.findMaxOverlap(stateOI, sol_vec[1].vec)
    transition_dipole_moment_vec = []
    mag = 0
    for (indF, field) in enumerate(x)
        stateC = Adiabatic ? (rev ? reverse(sol_vec)[indF].vec[:, startingInd] : sol_vec[indF].vec[:, startingInd])  : diabaticRamp(startingState, sol_vec, [x[1], field], Field = var)
        #println(typeof(stateC))
        #println(calculate.transition_dipole_moment(mol, stateC, rev ? reverse(sol_vec)[indF] : sol_vec[indF], comp = comp))
        tdm_field =  sum(abs.(calculate.transition_dipole_moment(mol, stateC, rev ? reverse(sol_vec)[indF].vec : sol_vec[indF].vec, comp = normalize!(comp))), dims = 2)
        if maximum(tdm_field) > mag
            mag = maximum(tdm_field)
        end

        push!(transition_dipole_moment_vec,tdm_field)
    end

    x = rev ? reverse(x) : x
    lplot = 0
    for rotationalStates in indsOI, state in (rotationalStates[1] + 1):rotationalStates[2]
        tdms = [tdm[state] for tdm in transition_dipole_moment_vec]
        color = rev ? log10.(reverse(tdms)) : log10.(tdms)
        #color = 1 .- color
    # Create an RGBA color map based on the normalized color values

        lplot = plotObj.lines!(ax, x, [sol_i.val[state] - energyRef for sol_i in sol_vec]*1e-9; color = color, colormap = plotObj.Reverse(:tokyo25), colorrange = (log10(1e-2), log10(1)) )

    end
    custom_formatter(values) = map(
        v -> "10" * plotObj.Makie.UnicodeFun.to_superscript(round(Int64, v)),
        values
    )
    
    cb = plotObj.Colorbar(fig[1, 2], lplot,   ticks =log10.([1e-2, 1e-1, 1e0]), tickformat=custom_formatter,  label = "Transition Dipole Moment", vertical = true)
    fig

end


end



end
