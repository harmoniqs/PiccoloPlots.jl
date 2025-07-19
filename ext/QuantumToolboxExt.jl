module QuantumToolboxExt

using PiccoloPlots
using QuantumToolbox
using Makie
using NamedTrajectories
using PiccoloQuantumObjects

using TestItems


function iso_to_bloch(ψ̃::AbstractVector{<:Real})
    s = QuantumObject(iso_to_ket(ψ̃))
    # Norm avoids Warning from QuantumToolbox
    return QuantumToolbox._state_to_bloch(s / norm(s))
end

function bloch_arrow(v, arrow_size)
    # Reduce size for correct arrow head position
    return (1 - arrow_size / norm(v)) * Vec3f(v...)
end

"""
    plot_bloch(
        traj::NamedTrajectory;
        state_name::Symbol = :ψ̃,
        kwargs...
    )

Plot the trajectory of a quantum state on the Bloch sphere.

TODO: Convert to density matrix in all cases.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to plot.

# Keyword Arguments
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.render`.

# Returns
- `fig`: The Makie `Figure` object.
"""
function QuantumToolbox.plot_bloch(
    traj::NamedTrajectory;
    index::Union{Nothing, Int} = nothing,
    state_name::Symbol = :ψ̃,
    kwargs...
)
    @assert state_name in traj.names
    bloch_vectors = iso_to_bloch.(eachcol(traj[state_name]))
    
    # Render Bloch sphere
    b = QuantumToolbox.Bloch()
    QuantumToolbox.add_points!(b, hcat(bloch_vectors...))

    fig, lscene = QuantumToolbox.render(b; kwargs...)

    if !isnothing(index)
        @assert 1 ≤ index ≤ length(bloch_vectors) "Invalid vector index."

        # Save for animation
        fig.attributes[:bloch] = b
        fig.attributes[:state_name] = state_name
        fig.attributes[:vec] = [bloch_arrow(bloch_vectors[index], b.vector_arrowsize[3])]

        # Draw the saved observable
        arrows!(
            lscene, [Point3f(0, 0, 0)], fig.attributes[:vec],
            linewidth = b.vector_width,
            arrowsize = Vec3f(b.vector_arrowsize...),
            rasterize = 3,
        )
    end

    return fig
end

function PiccoloPlots.plot_bloch!(
    fig::Figure,
    traj::NamedTrajectory,
    idx::Int;
    kwargs...
)
    @assert 1 ≤ idx ≤ traj.T "Invalid knot point index."

    if haskey(fig.attributes, :vec)
        b = fig.attributes[:bloch][]
        state_name = fig.attributes[:state_name][]
        v = iso_to_bloch(traj[idx][state_name])
        fig.attributes[:vec] = [bloch_arrow(v, b.vector_arrowsize[3])]
    end

    return fig
end

function PiccoloPlots.animate_bloch(
    traj::NamedTrajectory;
    fps::Int=30,
    mode::Symbol=:inline,
    filename="bloch_animation.mp4",
    kwargs...
)
    fig = plot_bloch(traj; index=1, kwargs...)

    return animate_figure(
        fig,
        1:traj.T,
        i -> plot_bloch!(fig, traj, i; kwargs...),
        mode=mode,
        fps=fps,
        filename=filename
    )
end


"""
    plot_wigner(
        traj::NamedTrajectory,
        idx::Int;
        state_name::Symbol=:ψ̃,
        xvec = -5:0.1:5,
        yvec = -5:0.1:5,
        projection::Val=Val(:two_dim),
        colorbar::Bool=true,
        kwargs...
    )

Plot the Wigner function of a quantum state in a trajectory.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `idx::Int`: The index of the state in the trajectory to plot.

# Keyword Arguments
- `name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `xvec`, `yvec`: Grids for plotting the Wigner function.

# Returns
- `fig`: The Makie `Figure`.
"""
function QuantumToolbox.plot_wigner(
    traj::NamedTrajectory,
    idx::Int;
    state_name::Symbol = :ψ̃,
    kwargs...
)
    @assert 1 ≤ idx ≤ traj.T "Invalid knot point index."
    state = QuantumObject(iso_to_ket(traj[idx][state_name]))
    fig, ax, hm = QuantumToolbox.plot_wigner(state; library = Val(:Makie), kwargs...)
    colsize!(fig.layout, 1, Aspect(1, 1.0))

    # Save attributes for animations
    fig.attributes[:state_name] = state_name
    fig.attributes[:ax] = ax
    fig.attributes[:hm] = hm
    fig.attributes[:label] = Label(fig[0, 1], "Timestep $idx", tellwidth=false)

    return fig
end

function PiccoloPlots.plot_wigner!(fig::Figure, traj::NamedTrajectory, idx::Int)
    @assert 1 ≤ idx ≤ traj.T "Invalid knot point index."

    # Extract attributes from the figure
    state_name = fig.attributes[:state_name][]
    hm = fig.attributes[:hm][]
    label = fig.attributes[:label][]

    state = QuantumObject(iso_to_ket(traj[idx][state_name]))
    W = transpose(wigner(state, hm[1][], hm[2][]))
    hm[3][] = W
    label.text[] = "Timestep $idx"
    return fig
end

function PiccoloPlots.animate_wigner(
    traj::NamedTrajectory;
    mode=:inline, 
    fps::Int=30,
    filename="wigner_animation.mp4",
    kwargs...
)
    # Setup: plot first frame and capture observables
    fig = plot_wigner(traj, 1; kwargs...)

    animate_figure(
        fig,
        1:traj.T,
        i -> plot_wigner!(fig, traj, i),
        mode=mode,
        fps=fps,
        filename=filename
    )

    return fig
end


#============================================================================#

@testitem "Test plot_bloch for Bloch sphere trajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    
    comps = (ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)), Δt = [1.0; 1.0],)
    traj = NamedTrajectory(comps)

    fig = plot_bloch(traj)

    @test fig isa Figure
end

@testitem "Test plot_bloch for Bloch sphere trajectory with one vector arrow shown" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    
    comps = (
        ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)),
        Δt = [1.0; 1.0],
    )

    traj = NamedTrajectory(comps)

    fig = plot_bloch(traj, index=1)

    @test fig isa Figure
end

@testitem "plot_bloch shows expected curved Bloch path" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    T = 20
    ts = range(0, π/2; length=T)

    kets = [
        QuantumObject(cos(θ) * [1.0 + 0im, 0.0 + 0im] + sin(θ) * [0.0 + 0im, 1.0 + 0im]) 
        for θ in ts
    ]

    iso_kets = ket_to_iso.(ψ.data for ψ in kets)

    ψ̃ = hcat(iso_kets...)
    Δt = fill(1.0, T)

    comps = (ψ̃ = ψ̃, Δt = Δt,)
    traj = NamedTrajectory(comps)
    fig = plot_bloch(traj, state_name=:ψ̃)
    @test fig isa Figure
end

@testitem "Plot Wigner function of coherent state" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloPlots
    using PiccoloQuantumObjects
    using CairoMakie

    N = 20
    α = 1.5 + 0.5im
    ψ = coherent(N, α)

    traj = NamedTrajectory((ψ̃ = hcat(ket_to_iso(ψ.data)), Δt = [1.0],))
    fig = plot_wigner(traj, 1, state_name=:ψ̃)

    @test fig isa Figure
end

end