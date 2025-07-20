module MakieAnimationsExt

using PiccoloPlots
using Makie
using NamedTrajectories


"""
    animate_figure(
        fig::Figure,
        frames::AbstractVector{Int},
        update_frame!::Function;
        mode::Symbol = :inline,
        fps::Int = 30,
        filename::String = "animation.mp4"
    )

Animate a Makie figure by updating frames.
"""
function PiccoloPlots.animate_figure(
    fig::Figure,
    frames::AbstractVector{Int},
    update_frame!::Function;
    mode::Symbol = :inline,
    fps::Int = 24,
    filename::String = "animation.mp4"
)
    if mode == :inline
        display(fig) # open the scene
        @async begin # don't block
            while isopen(fig.scene)
                for i in frames
                    update_frame!(i)
                    sleep(1 / fps)
                end
            end
        end
    elseif mode == :record
        record(fig, filename, frames; framerate=fps) do i
            update_frame!(i)
        end
    else
        error("Unsupported mode: $mode. Use :inline or :record.")
    end

    return fig
end


"""
    animate_name(
        traj::NamedTrajectory,
        name::Symbol;
        fps::Int = 24,
        mode::Symbol = :inline,
        fname::String = "name_animation.mp4",
        kwargs...
    )

Animate the evolution of a variable in a `NamedTrajectory`.
"""
function PiccoloPlots.animate_name(
    traj::NamedTrajectory,
    name::Symbol;
    fps::Int=24,
    mode::Symbol=:inline,
    filename="name_animation.mp4",
    kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set limits
    times = get_times(traj)
    xlims = [times[1], times[end]]
    ylims = collect(extrema(traj[name]))

    scatter!(ax, xlims, ylims, markersize=0.0)
    plot_name!(ax, traj, :a, indices=1:1, kwargs...)

    # TODO: Unclear how to set observables via indices, so redraw
    function update_frame!(i)
        empty!(ax.scene)
        scatter!(ax, xlims, ylims, markersize=0.0)
        plot_name!(ax, traj, :a, indices=1:i)
    end

    return animate_figure(
        fig,
        1:traj.T,
        update_frame!,
        mode=mode,
        fps=fps,
        filename=filename
    )
end


end