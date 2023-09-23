using Crystallography: AbstractCell, Lattice, edges, atomtypes, eachatom, eachatomgroup
using RecipesBase: @recipe, @series

@recipe function plot(lattice::Lattice; origin=zeros(eltype(lattice), 3))
    seriestype --> :path3d
    linewidth --> 1
    color --> :black
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    label := ""
    for edge in edges(lattice, origin)
        @series begin
            edge[:, 1], edge[:, 2], edge[:, 3]
        end
    end
end

@recipe function plot(
    cell::AbstractCell;
    origin=zeros(eltype(Lattice(cell)), 3),
    atomsconnected=false,
    edgewidth=1,
    bondwidth=2,
    bondcolor=:gray,
    atomsize=5,
)
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    lattice = Lattice(cell)
    @series begin
        origin := origin
        linewidth := edgewidth
        lattice
    end
    if atomsconnected
        allgroups = eachatomgroup(cell)
        for (i, groupᵢ) in enumerate(allgroups)
            for (j, groupⱼ) in enumerate(allgroups)
                if j > i
                    for (_, positionᵢₖ) in eachatom(groupᵢ)
                        for (_, positionⱼₗ) in eachatom(groupⱼ)
                            cpositionᵢₖ = lattice * positionᵢₖ  # Cartesian coordinates, do not use the same variable name as the loop variable!
                            cpositionⱼₗ = lattice * positionⱼₗ  # Cartesian coordinates, do not use the same variable name as the loop variable!
                            @series begin
                                seriestype := :path3d
                                linewidth := bondwidth
                                color --> bondcolor
                                label := ""
                                Tuple(map(collect, zip(cpositionᵢₖ, cpositionⱼₗ)))
                            end
                        end
                    end
                end
            end
        end
    end
    # Only show one label for each unique element
    for group in eachatomgroup(cell)
        coordinates = map(eachatom(group)) do (_, position)
            lattice * position  # Cartesian coordinates
        end
        seriestype --> :scatter3d
        markersize := atomsize
        markerstrokecolor --> :auto
        markerstrokewidth --> 0
        label := string(group.atom)
        XYZ = reduce(hcat, coordinates)
        @series begin
            XYZ[1, :], XYZ[2, :], XYZ[3, :]
        end
    end
end

@recipe function f(dispersion::DispersionRelation)
    points = interpolate(dispersion.paths)
    seriestype --> :path
    linewidth --> 1
    xguide --> ""
    yguide --> raw"energy"
    for band in bands
        @series begin
            reduce(hcat, band)
        end
    end
    for node in dispersion.path.nodes
        @series begin
            seriestype --> :vline
            seriescolor := :black  # Fix the axis color
            linewidth := 1  # This is an axis, don't change its width
            z_order --> :back
            label := ""
        end
    end
end
