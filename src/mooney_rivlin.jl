type Mooney_Rivlin
    C1::Float64
    C2::Float64
end


"""
Strain energy: W(C1, C2) = C1 * (J1 - 3) + C2 * (J2 - 3)

                 chain rule
stress: S = dWde    =         dWdJ1 * dJ1De + dWdJ2 * dJ2De + dWdJ3 * dJ3De

"""
function integrate_stress(model::Mooney_Rivlin, S, D, strain_voight, F)
    # The left Cauchy Green deformation_tensor
    # https://en.wikipedia.org/wiki/Finite_strain_theory#The_left_Cauchy.E2.80.93Green_or_finger_deformation_tensor
    C = F'*F
    J = det(F)

    # Invariants
    I1 = trace(C)
    I2 = 0.5 * (trance(C) ^ 2 - trace(C.^2))
    I3 = det(C)

    # Reduced invariants
    J1 = I1 * I3 ^ (-1/3.)
    J2 = I2 * I3 ^ (-2/3.)
    J3 = I3 ^ (1/2.)

    # Derivates of invariants
    # https://en.wikipedia.org/wiki/Tensor_derivative_%28continuum_mechanics%29#Derivatives_of_the_invariants_of_a_second-order_tensor
    dI1de = eye(3)
    dI2de = eye(3) - C'
    dI3de = det(C)*(inv(C))'

    dWdJ1 = C1
    dWdJ2 = C2
    dWdJ3 = K * (J3 - 1)

    dI1de = dI1De * I3^(-1/3) - 1/3 * I1*(I3^(-4/3))*dI3de
    dI2de = dI2De * I3^(-2/3) - 2/3 * I2*(I3^(-5/3))*dI3de
    dI3de = 1/2 * I3^(-1/2) * dI3de

    # Cauchy stress
    stress = 1 / J * (dWdJ1 * dJ1De + dWdJ2 * dJ2De + dWdJ3 * dJ3De) * F'
end
