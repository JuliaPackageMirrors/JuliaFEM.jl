using JuliaFEM


type Mooney_Rivlin
    C1::Float64
    C2::Float64
end


body = Problem(Elasticity, "block", 2)
body.properties.formulation = :plane_stress

nodes = Dict{Int64, Vector{Float64}}(
    1 => [0.0, 0.0],
    2 => [1.0, 0.0],
    3 => [1.0, 1.0],
    4 => [0.0, 1.0])

element = Element(Quad4, [1,2,3,4])
update!(element, "geometry", nodes)
element["youngs modulus"] = 9000.0
element["poissons ratio"] = 1/3
element["material model"] = Mooney_Rivlin(1.0, 1.0)

# dirichlet boundary conditions
boundary2 = Element(Seg2, [2,3])
update!(boundary2, "geometry", nodes)
boundary2["displacement traction force 2"] = [10.0, 0.0]

push!(body, element)
push!(body, boundary2)


# dirichlet boundary conditions
boundary_problem_1 = Problem(Dirichlet, "bc", 2, "displacement")
boundary = Element(Seg2, [1, 4])
update!(boundary, "geometry", nodes)
boundary["displacement 1"] = 0.0
boundary["displacement 2"] = 0.0
push!(boundary_problem_1, boundary)


solver = Solver("solve block problem")

push!(solver, body, boundary_problem_1)
call(solver)

disp = element("displacement", [1.0, 1.0], 0.0)
info("displacement at tip: $disp")
