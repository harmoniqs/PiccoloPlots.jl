module AnimatedPlots

export animate_figure

using Makie

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
function animate_figure(
    fig::Figure,
    frames::AbstractVector{Int},
    update_frame!::Function;
    mode::Symbol = :inline,
    fps::Int = 30,
    filename::String = "animation.mp4"
)
    if mode == :inline
        @async begin
            while isopen(fig)
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


end