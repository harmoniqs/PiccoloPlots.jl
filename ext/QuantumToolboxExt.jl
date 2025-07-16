module QuantumToolboxExt

using PiccoloPlots
using QuantumToolbox
using Makie
using NamedTrajectories
using PiccoloQuantumObjects
using TestItemRunner
using TestItems

"""
    plot_bloch(
        traj::NamedTrajectory;
        state_name::Symbol = :ψ̃,
        kwargs...
    )

Plot the trajectory of a quantum state on the Bloch sphere.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to plot.

# Keyword Arguments
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.render`.

# Returns
A tuple `(fig, lscene, states)` where:
- `fig`: The Makie `Figure` object.
- `lscene`: The 3D scene containing the Bloch sphere.
- `states`: The list of `QuantumObject`s plotted.
"""
function QuantumToolbox.plot_bloch(
    traj::NamedTrajectory;
    state_name::Symbol = :ψ̃,
    show_vector_at::AbstractVector{Int}=[],
    kwargs...
)
    iso_states = eachcol(traj[state_name])
    kets = iso_to_ket.(iso_states)
    states = QuantumObject.(kets)
    bloch_vectors = QuantumToolbox._state_to_bloch.(states)

    b = QuantumToolbox.Bloch()
    QuantumToolbox.add_points!(b, hcat(bloch_vectors...))

    for idx in show_vector_at
        @assert 1 ≤ idx ≤ traj.T "Invalid knot point index."
        QuantumToolbox.add_vectors!(b, bloch_vectors[idx])
    end

    fig, lscene = QuantumToolbox.render(b; kwargs...)

    display(fig)
    return fig, lscene, states
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
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `xvec`, `yvec`: Grids for plotting the Wigner function.
- `projection::Val`: Type of projection for visualization (`Val(:two_dim)` or `Val(:three_dim)`).
- `colorbar::Bool`: Whether to display a colorbar.
- `kwargs...`: Additional keyword arguments passed to the plot.

# Returns
A tuple `(fig, ax, hm, states)` where:
- `fig`: The Makie `Figure`.
- `ax`: The axis containing the plot.
- `hm`: The heatmap handle.
- `states`: The list of quantum states plotted.
"""
function QuantumToolbox.plot_wigner(
    traj::NamedTrajectory,
    idx::Int;
    state_name::Symbol = :ψ̃,
    xvec = -5:0.1:5,
    yvec = -5:0.1:5,
    projection::Val = Val(:two_dim),
    colorbar::Bool = true,
    kwargs...
)
    @assert 1 ≤ idx ≤ traj.T "Invalid knot point index."

    iso_states = eachcol(traj[state_name])
    kets = iso_to_ket.(iso_states)
    states = QuantumObject.(kets)

    fig, ax, hm = QuantumToolbox.plot_wigner(
        states[idx];
        library = Val(:Makie),
        xvec = xvec,
        yvec = yvec,
        projection = projection,
        colorbar = colorbar
    )

    return fig, ax, hm, states
end


#============================================================================#

@testitem "Test plot_bloch for Bloch sphere trajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    CairoMakie.activate!()
    x = ComplexF64[1.0; 0.0]
    y = ComplexF64[0.0, 1.0]
    
    comps = (
        ψ̃ = hcat(ket_to_iso(x), ket_to_iso(y)),
        Δt = [1.0; 1.0],
    )

    traj = NamedTrajectory(comps)

    fig, lscene, states = plot_bloch(traj)

    @test fig isa Figure
    @test lscene isa LScene
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

    fig, lscene, states = plot_bloch(traj, show_vector_at=[1])

    @test fig isa Figure
    @test lscene isa LScene
end

@testitem "plot_bloch shows expected curved Bloch path" begin
    using PiccoloPlots
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

    comps = (
        ψ̃ = ψ̃,
        Δt = Δt,
    )

    traj = NamedTrajectory(comps)

    fig, lscene = plot_bloch(traj, state_name=:ψ̃)

    @test isa(fig, Figure)
    @test isa(lscene, LScene)
end

@testitem "plot_bloch shows expected curved Bloch path with multiple vectors" begin
    using PiccoloPlots
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    T = 20
    ts = range(0, π/2; length=T)

    kets = [QuantumObject(cos(θ) * [1.0 + 0im, 0.0 + 0im] + sin(θ) * [0.0 + 0im, 1.0 + 0im]) for θ in ts]

    iso_kets = ket_to_iso.(ψ.data for ψ in kets)

    ψ̃ = hcat(iso_kets...)
    Δt = fill(1.0, T)

    comps = (
        ψ̃ = ψ̃,
        Δt = Δt,
    )

    traj = NamedTrajectory(comps)

    fig, lscene = PiccoloPlots.plot_bloch(Val(:Makie), traj, state_name=:ψ̃, show_vector_at=[1, 10, 15])

    @test isa(fig, Figure)
    @test isa(lscene, LScene)
end

@testitem "plot_bloch with problem-constructed trajectory" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects
    using QuantumCollocation
    using CairoMakie
    
    T = 50
    Δt = 0.2
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]])
    ψ_inits = Vector{ComplexF64}.([[1.0, 0.0], [0.0, 1.0]])
    ψ_targets = Vector{ComplexF64}.([[0.0, 1.0], [1.0, 0.0]])

    prob = QuantumStateSmoothPulseProblem(
        sys, ψ_inits, ψ_targets, T, Δt;
        piccolo_options=PiccoloOptions(verbose=false),
    )

    solve!(prob, max_iter=100, print_level=5)

    traj = prob.trajectory
    fig, lscene = plot_bloch(traj, state_name=:ψ̃1)

    @test fig isa Figure
    @test lscene isa LScene
end

@testitem "plot_bloch with problem-constructed trajectory and show_vector_at" begin
    using PiccoloPlots
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects
    using QuantumCollocation
    using CairoMakie
    
    T = 50
    Δt = 0.2
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]])
    ψ_inits = Vector{ComplexF64}.([[1.0, 0.0], [0.0, 1.0]])
    ψ_targets = Vector{ComplexF64}.([[0.0, 1.0], [1.0, 0.0]])

    prob = QuantumStateSmoothPulseProblem(
        sys, ψ_inits, ψ_targets, T, Δt;
    )

    solve!(prob, max_iter=100, print_level=5)

    traj = prob.trajectory
    fig, lscene = plot_bloch(traj, state_name=:ψ̃1, show_vector_at=[1, 10])

    @test fig isa Figure
    @test lscene isa LScene
end

@testitem "Plot Wigner function of coherent state" begin
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloPlots
    using PiccoloQuantumObjects
    using CairoMakie

    CairoMakie.activate!()

    N = 20
    α = 1.5 + 0.5im
    ψ = coherent(N, α)

    traj = NamedTrajectory((
        ψ̃ = hcat(ket_to_iso(ψ.data)), Δt = [1.0], 
    ))
    fig, ax, hm, _ = plot_wigner(traj, 1, state_name=:ψ̃)

    @test fig isa Figure
    @test ax isa Axis
    @test hm !== nothing
end

end