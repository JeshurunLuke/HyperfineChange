import Base.*
import Rotations: RotZYZ

struct MoleculeHamiltonian{T}
    MolOp
    Hhfs::T
    Hzeem::T
    Hdc::T
    Hac::Vector{T}
end
using WignerSymbols
using Base.Threads
using SparseArrays
using LinearAlgebra


function generateHamiltonian(Molecule::moleculeProperties; dirB= [0, 0, 1], dirE = [0, 0, 1.0], Beams = Dict())#, beams = OpticalTrap[])
    Hfs = Matrix(HyperFine_Ham(Molecule))
    Hzeem = Matrix(zeeman_ham(Molecule, normalize!(dirB)))
    Hdc = Matrix(dc(Molecule, normalize!(dirE)))
    Hac = Matrix{ComplexF64}[]#
    for (key, value) in Beams
        push!(Hac,  Matrix(AC_Stark(Molecule, value)))
    end
    if isempty(Hac)
        push!(Hac, zeros(ComplexF64, size(Hdc)))
    end
    MoleculeHamiltonian(Molecule, Hfs, Hzeem, Hdc, Hac)
end


function generateHamiltonianSparse(Molecule::moleculeProperties; dirB= [0, 0, 1], dirE = [0, 0, 1.0], Beams = Dict())#, beams = OpticalTrap[])
    Hfs = HyperFine_Ham(Molecule)
    Hzeem = zeeman_ham(Molecule, normalize!(dirB))
    Hdc = dc(Molecule, normalize!(dirE))
    Hac = SparseMatrixCSC{ComplexF64}[]#
    for (key, value) in Beams
        push!(Hac,  AC_Stark(Molecule, value))
    end
    if isempty(Hac)
        push!(Hac, spzeros(ComplexF64, size(Hdc)))
    end

    MoleculeHamiltonian(Molecule, Hfs, Hzeem, Hdc, Hac)
end


const PhysConstants = Dict("c" => 2.99798458E8, "hbar" => 1.054571E-34, "kb" => 1.3806452E-23, "e0" => 8.85418E-12, "q_c" => 1.60217663e-19, "m_e" => 9.11e-31)

function rotation_matrix_to_align_with_x_axis(v::Vector{T}) where T
    v_normalized = normalize(v)

    rotation_axis = cross(v_normalized, [1.0, 0.0, 0.0])

    if isapprox(norm(rotation_axis), 0.0)
        return Matrix{Float64}(I, 3, 3) # Return identity matrix
    end

    rotation_axis = normalize(rotation_axis)

    # Find the angle between the vector and the x-axis
    theta = acos(v_normalized[1])

    # Compute the rotation matrix using the Rodrigues' rotation formula
    K = [0 -rotation_axis[3] rotation_axis[2];
         rotation_axis[3] 0 -rotation_axis[1];
         -rotation_axis[2] rotation_axis[1] 0]

    rot_mat = I + sin(theta)*K + (1-cos(theta))*K^2

    return rot_mat
end


struct OpticalBeam
    kvec
    pol_k
    rotation_matrix
    pol_space
    Intensity
end





function compute_SHC_elements!(operator, q, N; order = 2)
    
    for (j, state) in enumerate(getBasis(N.spin)), (k, state2) in enumerate(getBasis(N.spin))
        if true #abs(state[2] - state2[2]) <= order && -1*state[1] + state2[1] + q == 0
            mN, N1 = state
            mNp, Np = state2 
            

            operator[j, k] = (-1)^(mN) * sqrt((2*N1 + 1)*(2*Np + 1)) * wigner3j(N1, order, Np, -mN, q, mNp) * wigner3j(N1, order, Np, 0, 0, 0)
        end
    end
end



function compute_quad_elements!(operator, q, node; order = 2)
    base = getBasis(node.spin)
    for (j, state) in enumerate(base), (k, state2) in enumerate(base)
        if true #abs(state[2] - state2[2]) <= order && -1*state[1] + state2[1] + q == 0
    
            mJ1, J1 = state 
            mJ2, J2 = state2
            operator[j, k] = (-1)^(J1-mJ1) * wigner3j(J1, order, J1, -mJ1, q, mJ2)/ wigner3j(J1, order, J1, -J1, 0, J1)
        end
    end
end

function quad_moment(mol::moleculeProperties)
    F = mol.basisTree
    EndNodes = endNode(F)
    nucleiQuad = []

    for (i, node) in enumerate(EndNodes[2:end])
        i += 1
        dimsA = i == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(i+ 1):end]])
        dimsB = i == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(i-1):-1:1]])

        Iquad = [spzeros(ComplexF64, length(node.spin), length(node.spin)) for _ in -2:2] 
        I_a = sparse(1.0I, dimsA, dimsA)
        I_b = sparse(1.0I, dimsB, dimsB)


        @threads for i = 1:5
            q = i - 3
            compute_quad_elements!(Iquad[i], q, node, order = 2)
        end
        
        push!(nucleiQuad, [kron(I_b, kron(sparse(ord), I_a)) for ord in Iquad])
    end
    nucleiQuad
end


function DipoleMatrix(mol::moleculeProperties)
    F = mol.basisTree
    inds = [node.name == "N" ? 1 : 0 for node in endNode(F)] 
    N = endNode(F)[inds .== 1][1]
    dims = prod([length(node.spin) for node in endNode(F)[inds .== 0]])
    dipoleOperator = [spzeros(ComplexF64, length(N.spin), length(N.spin)) for _ in -1:1] 
    @threads for i in 1:3
        q = i - 2
        compute_SHC_elements!(dipoleOperator[i], q, N, order = 1)
    end

    [kron(sparse(dip), sparse(1.0I, dims, dims)) for dip in dipoleOperator]
end




function getDipoleMatrix(mol::moleculeProperties)
    dipm, dip0, dipp =  DipoleMatrix(mol)
    dip = -1*mol.Constants["d0"]*[(dipm + adjoint(dipm))/sqrt(2), 1im*(dipm - adjoint(dipm))/sqrt(2), dip0]
end




function T2_C(mol::moleculeProperties)
    F = mol.basisTree
    inds = [node.name == "N" ? 1 : 0 for node in endNode(F)] 
    N = endNode(F)[inds .== 1][1]
    dims = prod([length(node.spin) for node in endNode(F)[inds .== 0]])
    T2C = [spzeros(ComplexF64, length(N.spin), length(N.spin)) for _ in -2:2] 
    @threads for i in 1:5
        q = i - 3
        compute_SHC_elements!(T2C[i], q, N, order = 2)
    end

    [kron(sparse(q), sparse(1.0I, dims, dims)) for q in T2C]
end




function makeT2(I1, I2)
    T2m2 = 0.5 * (I1[1] * I2[1] - 1im*I1[1] * I2[2] - 1im*I1[2] * I2[1] - I1[2] * I2[2])
    T2p2 = 0.5 * (I1[1] * I2[1] + 1im*I1[1] * I2[2] + 1im*I1[2] * I2[1] - I1[2] * I2[2])

    T2m1 = 0.5 * (I1[1] * I2[3] - 1im*I1[2] * I2[3] + I1[3] * I2[1] - 1im*I1[3] * I2[2])
    T2p1 = -0.5 * (I1[1] * I2[3] + 1im*I1[2] * I2[3] + I1[3] * I2[1] + 1im*I1[3] * I2[2])

    T20 = -sqrt(1/6) * (I1[1] * I2[1] + I1[2] * I2[2]) + sqrt(2/3) * I1[3] * I2[3]

    T = [T2m2, T2m1, T20, T2p1, T2p2]

    return T
end


function tensor_nuclear(M::moleculeProperties)
    T1 = T2_C(M)
    T2 = makeT2(M.I1, M.I2)
    sqrt(6)*M.Constants["C3"]*tensor_dot(T1, T2, order = 2)
end




function dc(M::moleculeProperties, dirE)
    sum(dirE.* getDipoleMatrix(M))
end



function AC_Stark(M::moleculeProperties, Beam::OpticalBeam)
    h = PlanckConstant.val
    println(Beam.Intensity)
    mag = sqrt(2*Beam.Intensity/(PhysConstants["c"]*PhysConstants["e0"]))
    a0, a2 = M.Constants["a0"], M.Constants["a2"]
    dims = prod([length(node.spin) for node in endNode(M.basisTree)])
    pol = normalize(Beam.pol_space)
    
    A2 = sqrt(3/2)*T2_C(M)
    P2 = [scaler*sparse(1.0I, dims, dims) for scaler in makeT2(pol, conj.(pol))]

    -mag^2*(sparse(1.0I, dims, dims)*a0 + tensor_dot(A2, P2, order = 2)*a2)/4
end


OpticalBeam(kvec, pol_k) = OpticalBeam(kvec, pol_k, 1.0)
function OpticalBeam(kvec, pol_k, Intensity)
    kvec, pol_k = normalize(kvec), normalize(pol_k)
    
    rot_mat = rotation_matrix_to_align_with_x_axis(kvec)
    pol_space = adjoint(rot_mat)*pol_k ## Reasonable>?
    if !isapprox([1.0, 0, 0], rot_mat*kvec)
        throw("Improper Rotation: $(rot_mat*kvec) instead of [1.0, 0.0, 0.0]")
    end
    
    if pol_k[1] != 0.0
        throw("$pol_k Improper Polarization: (pol_k should not have x components)")
    end


    if  !isapprox(dot(kvec, pol_space), 0.0, atol=1e-8)
        throw("$pol_space Kvec is not perpendicular to the polarization")
    end


    OpticalBeam(kvec, pol_k, rot_mat, pol_space, Intensity)
end






*(a::Vector{T}, b::Vector{T}) where {T<:AbstractSparseMatrix} = sum([a[i]*b[i] for i = 1:size(a)[1]])
function generateRotational(Nrot, Brot, Drot)
    Nsq = Nrot*Nrot
    Brot*Nsq - Drot*Nsq*Nsq
end
function HyperFine_Ham(M::moleculeProperties)
    consts = M.Constants
    scalarInteractions = generateRotational(M.N[1:3], consts["Brot"] , consts["Drot"]) + scalar_nuclear(consts["C1"], M.N[1:3], M.I1[1:3]) + scalar_nuclear(consts["C2"], M.N[1:3], M.I2[1:3]) + scalar_nuclear(consts["C4"], M.I1[1:3], M.I2[1:3])
    tensorInteraction =  quadrapole(M)  #+ tensor_nuclear(M) 
    
    
    scalarInteractions + tensorInteraction 
end

function quadrapole(mol::moleculeProperties)
    Q1, Q2 = mol.Constants["Q1"], mol.Constants["Q2"]
    TdE = T2_C(mol)
    Tq2, Tq1 = quad_moment(mol)
    Hq = (Q1*tensor_dot(Tq1, TdE, order = 2) + Q2*tensor_dot(Tq2, TdE, order = 2))/4

end

function tensor_dot(T1, T2; order = 2)
    Tprod = spzeros(ComplexF64, size(T1[1])...)
    T2 = reverse(T2)
    for (i, q) in enumerate(-order:order)
        Tprod += ((-1)^q)*T1[i]*T2[i]
    end
    Tprod
end




function zeeman_ham(M::moleculeProperties, Bvec)
    consts = M.Constants
    Nproj, I1proj, I2proj = sum(Bvec.*M.N[1:3]), sum(Bvec.*M.I1[1:3]), sum(Bvec.*M.I2[1:3])
    H = -consts["Mu1"]*I1proj -consts["Mu2"]*I2proj -consts["MuN"]*Nproj
end


function scalar_nuclear(consts, J1, J2)
    consts*J1*J2
end


