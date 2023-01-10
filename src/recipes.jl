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
