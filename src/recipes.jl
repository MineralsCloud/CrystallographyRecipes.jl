using RecipesBase: @recipe, @series

@recipe function plot(lattice::Lattice)
    seriestype --> :path3d
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    color --> :grey
    tick_direction --> :out
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    legend --> :none
    frame --> :box
    for edge in edges(lattice)
        @series begin
            edge[:, 1], edge[:, 2], edge[:, 3]
        end
    end
end

@recipe function plot(cell::Cell)
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    tick_direction --> :out
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    frame --> :box
    for edge in edges(cell.lattice)
        @series begin
            seriestype --> :path3d
            color --> :grey
            label --> :none
            edge[:, 1], edge[:, 2], edge[:, 3]
        end
    end
    transform = CartesianFromFractional(cell.lattice)
    for (atom, position) in zip(cell.atoms, cell.positions)
        position = transform(position)
        @series begin
            seriestype --> :scatter3d
            markerstrokecolor --> :auto
            markerstrokewidth --> 0
            label --> string(atom)
            Tuple(Base.vect(coordinate) for coordinate in position)
        end
    end
end
