module CrystallographyRecipes

using CrystallographyBase: AbstractLattice
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

@recipe function plot(lattice::AbstractLattice)
    seriestype --> :path3d
    linewidth --> 1
    color --> :black
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    label := ""
    for edge in edges(lattice)
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

@userplot BandsPlot
@recipe function f(plot::BandsPlot; specialpoints=[[]])
    if !(eltype(specialpoints) <: AbstractVector)
        throw(ArgumentError("`specialpoints` must be a vector of vectors!"))
    end
    dispersions, recip_lattice = plot.args
    paths = collect(dispersion.path for dispersion in dispersions)
    bands = reduce(vcat, (dispersion.bands for dispersion in dispersions))
    xticks = cumsum(normalize_lengths(paths, recip_lattice))
    prepend!(xticks, zero(eltype(xticks)))  # Don't forget the initial point!
    xticklabels = _makexticklabels(specialpoints)
    xticks --> (xticks, xticklabels)
    xlims --> extrema(xticks)
    xguide --> "k"
    for band in eachcol(bands)
        @series begin
            seriestype --> :path
            linewidth --> 1
            label --> ""
            eachindex(band) ./ length(band), band
        end
    end
    split = cumsum(length.(specialpoints))[begin:(end - 1)]  # Do not include the last point
    for (i, xtick) in enumerate(xticks)
        seriestype --> :vline
        seriescolor := :black  # Fix the axis color
        linewidth := 1  # This is an axis, don't change its width
        z_order --> :back
        label := ""
        if i in split
            @series begin
                linestyle --> :solid
                [xtick]
            end
        else
            @series begin
                linestyle --> :dot
                [xtick]
            end
        end
    end
end

function _makexticklabels(specialpoints)
    return collect(
        Iterators.flatten(
            map(firstindex(specialpoints):lastindex(specialpoints)) do i
                points = specialpoints[i]
                labels = String[]
                for (j, point) in enumerate(points)
                    if j == firstindex(points)
                        if i == firstindex(specialpoints)
                            push!(labels, string(point))
                        else
                            continue
                        end
                    elseif j == lastindex(points)
                        if i == lastindex(specialpoints)
                            push!(labels, string(point))
                        else
                            push!(labels, string(point, '|', specialpoints[i + 1][begin]))
                        end
                    else
                        push!(labels, string(point))
                    end
                end
                labels
            end,
        ),
    )
end

@userplot DOSPlot
@recipe function f(plot::DOSPlot; vertical=false)
    energies, dos = plot.args
    I = sortperm(energies)  # Remember to sort!
    yguide --> "DOS"
    seriestype --> :path
    xlims --> extrema(energies)
    ylims --> extrema(dos)
    if vertical
        permute --> (:x, :y)  # See https://discourse.julialang.org/t/is-there-a-simple-way-to-swap-the-x-and-y-axes-in-plots-jl-after-a-figure-is-plotted/103599/5
    end
    return energies[I], dos[I]
end

end
