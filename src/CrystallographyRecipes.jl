module CrystallographyRecipes

using CrystallographyBase: Lattice, Cell, CartesianFromFractional, latticevectors

export vertices, edge, edges

function vertices(lattice::Lattice)
    O = zeros(eltype(lattice), 3)
    A, B, C = latticevectors(lattice)
    A_B, A_C, B_C, ABC = A + B, A + C, B + C, A + B + C
    return [O, A, B, C, A_B, A_C, B_C, ABC]
end

edge(A, B) = hcat(([Aᵢ, Bᵢ] for (Aᵢ, Bᵢ) in zip(A, B))...)

function edges(lattice::Lattice)
    O, A, B, C, AB, AC, BC, ABC = vertices(lattice)
    return [
        edge(O, A),
        edge(O, C),
        edge(A, AC),
        edge(C, AC),
        edge(O, B),
        edge(A, AB),
        edge(C, BC),
        edge(AC, ABC),
        edge(B, AB),
        edge(B, BC),
        edge(AB, ABC),
        edge(BC, ABC),
    ]
end

include("recipes.jl")

end
