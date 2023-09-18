using CrystallographyBase: Lattice, Cell, CartesianFromFractional, edges
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
    lattice = Lattice(cell)
    @series begin
        label --> :none
        lattice
    end
    transform = CartesianFromFractional(cell.lattice)
    # Only show one label for each unique element
    types = unique(string.(cell.atoms))
    for type in types
        count = 0
        indices = findall(atom -> string(atom) == type, cell.atoms)
        for position in cell.positions[indices]
            position = transform(position)
            @series begin
                seriestype --> :scatter3d
                markerstrokecolor --> :auto
                markerstrokewidth --> 0
                label := count >= 1 ? "" : type
                Tuple(Base.vect(coordinate) for coordinate in position)
            end
            count += 1
        end
    end
end
