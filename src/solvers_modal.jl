# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/JuliaFEM.jl/blob/master/LICENSE.md

"""
Examples
--------

julia> problems = get_problems()
julia> solver = Solver(Modal)
julia> push!(solver, problems...)
julia> solver()

"""
type Modal <: AbstractSolver
    geometric_stiffness :: Bool
    eigvals :: Vector
    eigvecs :: Matrix
    nev :: Int
    which :: Symbol
end

function Modal(nev=10, which=:SM)
    solver = Modal(false, Vector(), Matrix(), nev, which)
end

function call(solver::Solver{Modal}; show_info=true, debug=false, bc_invertible=false)
    show_info && info(repeat("-", 80))
    show_info && info("Starting natural frequency solver")
    show_info && info("Increment time t=$(round(solver.time, 3))")
    show_info && info(repeat("-", 80))
    initialize!(solver)
    # assemble all field problems
    info("Assembling problems ...")
    tic()
    for problem in get_field_problems(solver)
        assemble!(problem, solver.time)
        assemble!(problem, solver.time, Val{:mass_matrix})
    end
    for problem in get_boundary_problems(solver)
        assemble!(problem, solver.time)
    end
    t1 = round(toq(), 2)
    info("Assembled in $t1 seconds.")
    M, K, Kg, f = get_field_assembly(solver)
    Kb, C1, C2, D, fb, g = get_boundary_assembly(solver)
    K = K + Kb
    f = f + fb
    if solver.properties.geometric_stiffness
        K += Kg
    end

    @assert nnz(D) == 0
    @assert C1 == C2

    tic()

    if bc_invertible
        P, h = create_projection(C1, g, Val{:invertible})
    else
        P, h = create_projection(C1, g)
    end
    K_red = P'*K*P
    M_red = P'*M*P
    # make sure matrices are symmetric
    K_red = 1/2*(K_red + K_red')
    M_red = 1/2*(M_red + M_red')

#=
    ndim = size(C1,1)
    nz = get_nonzero_rows(C1)
    nz = setdiff(collect(1:ndim), nz)
    g = zeros(ndim)
    P = spzeros(ndim, ndim)
    for j in nz
        P[j,j] = 1.0
    end
    K_red = P'*K*P
    M_red = P'*M*P
    # make sure matrices are symmetric
    K_red = 1/2*(K_red + K_red')
    M_red = 1/2*(M_red + M_red')

    #=
    K_red = K[nz,nz]
    M_red = M[nz,nz]
    # make sure matrices are symmetric
    K_red = 1/2*(K_red + K_red')
    M_red = 1/2*(M_red + M_red')
    =#

    #=
    K_red = copy(K)
    M_red = copy(M)
    for j=1:size(K_red)
        j in nz && continue
        K_red[j,:] = 0.0
        K_red[:,j] = 0.0
        M_red[j,:] = 0.0
        M_red[:,j] = 0.0
    end
    =#
=#

    t1 = round(toq(), 2)
    info("Eliminated dirichlet boundaries in $t1 seconds.")

    nz = get_nonzero_rows(K_red)
    ndofs = solver.ndofs
    props = solver.properties
    info("Calculate $(props.nev) eigenvalues...")
    if debug && length(nz) < 100
        info("Stiffness matrix:")
        dump(round(full(K[nz, nz])))
        info("Mass matrix:")
        dump(round(full(M[nz, nz])))
    end

    tic()
    om2 = nothing
    X = nothing
    try
        om2, X = eigs(K_red[nz,nz], M_red[nz,nz]; nev=props.nev, which=props.which)
    catch
        info("failed to calculate eigenvalues")
        info("reduced system")
        info("is K symmetric? ", issym(K_red[nz,nz]))
        info("is M symmetric? ", issym(M_red[nz,nz]))
        info("is K positive definite? ", isposdef(K_red[nz,nz]))
        info("is M positive definite? ", isposdef(M_red[nz,nz]))
        k1 = maximum(abs(K_red[nz,nz] - K_red[nz,nz]'))
        m1 = maximum(abs(M_red[nz,nz] - M_red[nz,nz]'))
        info("K 'skewness' (max(abs(K - K'))) = ", k1)
        info("M 'skewness' (max(abs(M - M'))) = ", m1)

        info("original matrix")
        info("is K symmetric? ", issym(K[nz,nz]))
        info("is M symmetric? ", issym(M[nz,nz]))
        info("is K positive definite? ", isposdef(K[nz,nz]))
        info("is M positive definite? ", isposdef(M[nz,nz]))
        k1 = maximum(abs(K[nz,nz] - K[nz,nz]'))
        m1 = maximum(abs(M[nz,nz] - M[nz,nz]'))
        info("K 'skewness' (max(abs(K - K'))) = ", k1)
        info("M 'skewness' (max(abs(M - M'))) = ", m1)

        rethrow()
    end
    info("Eigenvalues computed in $t1 seconds. Eigenvalues: $om2")

    props.eigvals = om2
    props.eigvecs = zeros(ndofs, length(om2))
    v = zeros(ndofs)
    for i=1:length(om2)
        fill!(v, 0.0)
        v[nz] = X[:,i]
        props.eigvecs[:,i] = P*v + g
    end
    t1 = round(toq(), 2)

    for i=1:length(om2)
        freq = real(sqrt(om2[i])/(2.0*pi))
        u = props.eigvecs[:,i]
        field_dim = get_unknown_field_dimension(solver)
        field_name = get_unknown_field_name(solver)
        if field_dim != 1
            nnodes = round(Int, length(u)/field_dim)
            u = reshape(u, field_dim, nnodes)
            u = Vector{Float64}[u[:,i] for i in 1:nnodes]
        end
        for problem in get_problems(solver)
            for element in get_elements(problem)
                connectivity = get_connectivity(element)
                update!(element, field_name, freq => u[connectivity])
            end
        end
    end

    return true
end

