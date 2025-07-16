module GLMakieExt

using PiccoloPlots
using QuantumToolbox
using NamedTrajectories
using PiccoloQuantumObjects
using TestItemRunner
using TestItems
using GLMakie
using Colors

function PiccoloPlots.animate_bloch(
    traj::NamedTrajectory;
    fps::Int=30,
    repeat::Bool=false,
    kwargs...
)

    fig, lscene, states = plot_bloch(traj; kwargs...)
    bloch_vectors = QuantumToolbox._state_to_bloch.(states)

    bx = getindex.(bloch_vectors, 1)
    by = getindex.(bloch_vectors, 2)
    bz = getindex.(bloch_vectors, 3)

    index = Observable(1)
    moving_point = @lift(Point3f(bx[$index], by[$index], bz[$index]))
    scatter!(lscene, moving_point, color=:red, markersize=15)

    display(fig)

    @async begin
        disp = true
        while (!(fig.scene.isclosed)) && disp
            for i in 1:length(bloch_vectors)
                if fig.scene.isclosed
                    break
                end

                index[] = i
                sleep(1 / fps)

                disp = disp && repeat
            end
        end
    end

    return fig, lscene, states
end


function PiccoloPlots.animate_name(
    traj::NamedTrajectory,
    name::Symbol;
    fps::Int=30,
    repeat::Bool=false,
)
    
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_name!(ax, traj, name)

    index = Observable(1)
    times = get_times(traj)
    data = traj[name]
    for i in 1:size(data, 1)
        moving_dot = @lift(Point2f(times[$index], data[i, $index]))
        scatter!(ax, moving_dot, color=:red, markersize=10)
    end

    display(fig)

    @async begin
        disp = true
        while (!(fig.scene.isclosed)) && disp
            for t in 1:traj.T
                if fig.scene.isclosed
                    break
                end

                index[] = t
                sleep(1 / fps)

                disp = disp && repeat
            end
        end
    end

    return fig, ax
end

function PiccoloPlots.animate_wigner(
    traj::NamedTrajectory;
    fps::Int = 30,
    repeat::Bool=false,
    kwargs...
)

    fig, ax, hm, states = PiccoloPlots.plot_wigner(traj; kwargs...)
    
    label = Label(fig[0, 1], "Timestep 1", tellwidth=false)
    
    display(fig)

    @async begin
        disp = true
        while (!(fig.scene.isclosed)) && disp
            for (i, ψ) in enumerate(states)
                if fig.scene.isclosed
                    break
                end

                label.text[] = "Timestep $i"
                W = transpose(wigner(ψ, xvec, yvec))
                hm[3][] = W 
                sleep(1 / fps)

                disp = disp && repeat
            end
        end
    end

    return fig, ax, states
end

end