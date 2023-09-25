using Crystallography:
    AbstractCell,
    Lattice,
    ShiftedLattice,
    DispersionRelation,
    edges,
    atomtypes,
    eachatom,
    eachatomgroup,
    normalize_lengths
using RecipesBase: @recipe, @userplot, @series

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
@recipe function plot(lattice::ShiftedLattice)
    seriestype --> :path3d
    linewidth --> 1
    color --> :black
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    label := ""
    for edge in edges(Lattice(lattice.original .+ lattice.by), lattice.by)
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

@userplot DispersionRelationsPlot
@recipe function f(plot::DispersionRelationsPlot; specialpoints=[], split=[])
    xguide --> "k"
    dispersions, recip_lattice = plot.args
    paths = collect(dispersion.path for dispersion in dispersions)
    bands = reduce(vcat, (dispersion.bands for dispersion in dispersions))
    𝐋 = accumulate(+, normalize_lengths(paths, recip_lattice))
    xlims --> extrema(𝐋)
    xticks --> (𝐋, string.(specialpoints))
    for band in eachcol(bands)
        @series begin
            seriestype --> :path
            linewidth --> 1
            label --> ""
            eachindex(band) ./ length(band), band
        end
    end
    for (i, L) in enumerate(𝐋)
        if i in split
            @series begin
                seriestype --> :vline
                seriescolor := :black  # Fix the axis color
                linewidth := 1  # This is an axis, don't change its width
                z_order --> :back
                label := ""
                [L]
            end
        end
    end
end
