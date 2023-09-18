using Crystallography:
    AbstractCell, Lattice, edges, atomtypes, eachatom, eachatomgroup, supercell
using RecipesBase: @recipe, @series

@recipe function plot(lattice::Lattice)
    seriestype --> :path3d
    linewidth --> 1
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    color --> :black
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    label --> :none
    for edge in edges(lattice)
        @series begin
            edge[:, 1], edge[:, 2], edge[:, 3]
        end
    end
end

@recipe function plot(cell::AbstractCell)
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    lattice = Lattice(cell)
    @series begin
        lattice
    end
    # Only show one label for each unique element
    for group in eachatomgroup(cell)
        coordinates = map(eachatom(group)) do (_, position)
            lattice * position  # Cartesian coordinates
        end
        @series begin
            seriestype --> :scatter3d
            markersize --> 5
            markerstrokecolor --> :auto
            markerstrokewidth --> 0
            label --> string(group.atom)
            map(Base.Fix2(getindex, 1), coordinates),
            map(Base.Fix2(getindex, 2), coordinates),
            map(Base.Fix2(getindex, 3), coordinates)
        end
    end
end
