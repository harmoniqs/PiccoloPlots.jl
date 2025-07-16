module GLMakieExt

using PiccoloPlots
using QuantumToolbox
using Makie
using NamedTrajectories
using PiccoloQuantumObjects
using TestItemRunner
using TestItems
using GLMakie
using Colors


function PiccoloPlots.animate_bloch(
    traj::NamedTrajectory;
    state_name::Symbol=:ψ̃,
    fps::Int=30,
    show_vector_at::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    repeat::Bool=false,
)

    fig, lscene, states = PiccoloPlots.plot_bloch(
        Val(:Makie), traj; state_name=state_name, show_vector_at=show_vector_at
    )
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
    traj::NamedTrajectory;
    state_name::Symbol=:x,
    fps::Int=30,
    repeat::Bool=false,
)
    
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_name!(ax, traj, state_name)

    index = Observable(1)
    times = get_times(traj)
    data = traj[state_name]
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
    state_name::Symbol = :ψ̃,
    fps::Int = 30,
    xvec = -5:0.1:5,
    yvec = -5:0.1:5,
    projection::Val = Val(:two_dim),
    colorbar::Bool = true,
    repeat::Bool=false,
)

    fig, ax, hm, states = PiccoloPlots.plot_wigner(
        Val(:Makie),
        traj;
        state_name = state_name,
        xvec = xvec,
        yvec = yvec,
        projection = projection,
        colorbar = colorbar
    )
    
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


# @testitem "Animate Bloch sphere with quarter-turn trajectory" begin
#     using PiccoloPlots
#     using QuantumToolbox
#     using NamedTrajectories
#     using PiccoloQuantumObjects
#     using GLMakie

#     GLMakie.activate!()

#     T = 20
#     ts = range(0, π/2; length=T)
#     kets = [QuantumObject(cos(θ)*[1.0, 0.0] + sin(θ)*[0.0, 1.0]) for θ in ts]
#     iso_kets = ket_to_iso.(ψ.data for ψ in kets)

#     traj = NamedTrajectory((
#         ψ̃ = hcat(iso_kets...),
#         Δt = fill(1.0, T),
#     ))

#     fig, lscene, states = PiccoloPlots.animate_bloch(traj; fps=10)

#     @test fig isa Figure
#     @test lscene isa LScene
#     @test length(states) == T
# end

# @testitem "Animate Bloch sphere with quarter-turn trajectory showing vectors" begin
#     using PiccoloPlots
#     using QuantumToolbox
#     using NamedTrajectories
#     using PiccoloQuantumObjects
#     using GLMakie

#     GLMakie.activate!()

#     T = 20
#     ts = range(0, π/2; length=T)
#     kets = [QuantumObject(cos(θ)*[1.0, 0.0] + sin(θ)*[0.0, 1.0]) for θ in ts]
#     iso_kets = ket_to_iso.(ψ.data for ψ in kets)

#     traj = NamedTrajectory((
#         ψ̃ = hcat(iso_kets...),
#         Δt = fill(1.0, T),
#     ))

#     fig, lscene, states = PiccoloPlots.animate_bloch(traj; fps=10, show_vector_at=[1, 10, 20])

#     @test fig isa Figure
#     @test lscene isa LScene
#     @test length(states) == T
# end


# @testitem "Animate scalar data with quarter-turn trajectory" begin
#     using PiccoloPlots
#     using NamedTrajectories
#     using GLMakie

#     GLMakie.activate!()

#     T = 20
#     ts = range(0, 2π; length=T)
#     sin_values = sin.(ts)
#     traj = NamedTrajectory((
#         x = reshape(sin_values, 1, :),
#         Δt = fill(1.0, T),
#     ))

#     fig, ax = PiccoloPlots.animate_name(traj, state_name=:x, fps=10)

#     @test fig isa Figure
#     @test ax isa Axis
# end


# @testitem "Animate Wigner function with quarter-turn trajectory" begin
#     using PiccoloPlots
#     using QuantumToolbox
#     using NamedTrajectories
#     using PiccoloQuantumObjects: ket_to_iso
#     using GLMakie

#     GLMakie.activate!()

#     T = 20
#     ts = range(0, π/2; length=T)
#     kets = [QuantumObject(cos(θ)*[1.0, 0.0] + sin(θ)*[0.0, 1.0]) for θ in ts]
#     iso_kets = ket_to_iso.(ψ.data for ψ in kets)

#     traj = NamedTrajectory((
#         ψ̃ = hcat(iso_kets...),
#         Δt = fill(1.0, T),
#     ))

#     fig, ax, states = PiccoloPlots.animate_wigner(traj; fps=5)

#     @test fig isa Figure
#     @test ax isa Axis
#     @test length(states) == T
# end

end