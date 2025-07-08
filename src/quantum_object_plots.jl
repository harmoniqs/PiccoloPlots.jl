module QuantumObjectPlots

export plot_unitary_populations
export plot_bloch
export plot_first_bloch
export plot_last_bloch
export plot_bloch_traj
export animate_bloch_old
export animate_bloch
export animate_name
export plot_wigner
export animate_wigner

using LaTeXStrings
# only need MakieCore for plot name
using MakieCore
using LinearAlgebra
using GeometryBasics
using NamedTrajectories
using PiccoloQuantumObjects
using TestItemRunner
using TestItems

function get_layout(index::Int, n::Int)
    √n = isqrt(n) + 1
    return ((index - 1) ÷ √n + 1, (index - 1) % √n + 1)
end

"""
    plot_unitary_populations(
        traj::NamedTrajectory;
        unitary_columns::AbstractVector{Int}=1:2,
        unitary_name::Symbol=:Ũ⃗,
        control_name::Symbol=:a,
        kwargs...
    )

Plot the populations of the unitary columns of the unitary matrix in the trajectory. `kwargs` are passed to [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/).

# Keyword Arguments
- `unitary_columns::AbstractVector{Int}`: The columns of the unitary matrix to plot the populations of. Default is `1:2`.
- `unitary_name::Symbol`: The name of the unitary matrix in the trajectory. Default is `:Ũ⃗`.
- `control_name::Symbol`: The name of the control in the trajectory. Default is `:a`.
- `kwargs...`: Additional keyword arguments passed to [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/).
"""
function plot_unitary_populations(
    traj::NamedTrajectory;
    unitary_columns::AbstractVector{Int}=1:2,
    unitary_name::Symbol=:Ũ⃗,
    control_name::Symbol=:a,
    kwargs...
)

    transformations = [
        (unitary_name => x -> abs2.(iso_vec_to_operator(x)[:, i]))
        for i ∈ unitary_columns
    ]

    transformation_labels = [
        L"P" for i ∈ unitary_columns
    ]

    transformation_titles = [
        L"Populations: $\left| U_{:, %$(i)}(t) \right|^2$" for i ∈ unitary_columns
    ]
    
    MakieCore.plot(traj, [control_name];
        transformations=transformations,
        transformation_titles=transformation_titles,
        transformation_labels=transformation_labels,
        kwargs...
    )
end

"""
    plot_bloch_traj(
        ::Val{:Makie},
        traj::NamedTrajectory;
        state_name::Symbol = :ψ̃,
        kwargs...
    )

Plot the trajectory of a quantum state on the Bloch sphere.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states to plot.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.

# Keyword Arguments
- `kwargs...`: Additional keyword arguments passed to `QuantumToolbox.render`.

# Returns
A tuple `(fig, lscene, states)` where:
- `fig`: The Makie `Figure` object.
- `lscene`: The 3D scene containing the Bloch sphere.
- `states`: The list of `QuantumObject`s plotted.
"""
function plot_bloch_traj end

"""
    plot_first_bloch(traj::NamedTrajectory; state_name::Symbol=:ψ̃, kwargs...)

Plot the first quantum state in a trajectory on the Bloch sphere.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.

# Returns
A tuple `(fig, lscene)` where:
- `fig`: The Makie `Figure` object.
- `lscene`: The 3D scene containing the Bloch sphere.
"""
function plot_first_bloch end

"""
    plot_last_bloch(traj::NamedTrajectory; state_name::Symbol=:ψ̃, kwargs...)

Plot the last quantum state in a trajectory on the Bloch sphere.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.

# Returns
A tuple `(fig, lscene)` where:
- `fig`: The Makie `Figure` object.
- `lscene`: The 3D scene containing the Bloch sphere.
"""
function plot_last_bloch end

"""
    plot_wigner(
        traj::NamedTrajectory;
        state_name::Symbol=:ψ̃,
        library::Val = Val(:Makie),
        xvec = -5:0.1:5,
        yvec = -5:0.1:5,
        projection::Val=Val(:two_dim),
        colorbar::Bool=true,
        kwargs...
    )

Plot the Wigner function of a quantum state in a trajectory.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `library::Val`: The plotting library to use (default `Val(:Makie)`).

# Keyword Arguments
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
function plot_wigner end

"""
    animate_bloch_old(traj::NamedTrajectory; state_name::Symbol=:ψ̃, fps::Int=60)

Animate the evolution of a quantum trajectory on the Bloch sphere without using Quantum Toolbox.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `fps::Int`: Frames per second for the animation. Default is `60`.

# Returns
The Makie `Figure` object displaying the animation.
"""
function animate_bloch_old end

"""
    animate_bloch(traj::NamedTrajectory; state_name::Symbol=:ψ̃, fps::Int=30)

Animate the evolution of a quantum trajectory on the Bloch sphere.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `fps::Int`: Frames per second for the animation. Default is `30`.

# Returns
A tuple `(fig, lscene, states)` where:
- `fig`: The Makie `Figure` object.
- `lscene`: The 3D scene containing the Bloch sphere.
- `states`: The list of quantum states.
"""
function animate_bloch end

"""
    animate_name(traj::NamedTrajectory; state_name::Symbol=:x, fps::Int=30)

Animate the evolution of a scalar or vector-valued variable in the trajectory.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing the variable.
- `state_name::Symbol`: The name of the variable in the trajectory. Default is `:x`.
- `fps::Int`: Frames per second for the animation. Default is `30`.

# Returns
A tuple `(fig, ax)` where:
- `fig`: The Makie `Figure`.
- `ax`: The axis containing the animated plot.
"""
function animate_name end

"""
    animate_wigner(
        traj::NamedTrajectory;
        state_name::Symbol=:ψ̃,
        fps::Int=30,
        xvec = -5:0.1:5,
        yvec = -5:0.1:5,
        projection::Val=Val(:two_dim),
        colorbar::Bool=true,
        kwargs...
    )

Animate the evolution of the Wigner function for a quantum trajectory.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `fps::Int`: Frames per second for the animation. Default is `30`.

# Keyword Arguments
- `xvec`, `yvec`: Grids for the Wigner function.
- `projection::Val`: Type of projection for visualization (`Val(:two_dim)` or `Val(:three_dim)`).
- `colorbar::Bool`: Whether to display a colorbar.
- `kwargs...`: Additional keyword arguments passed to the plot.

# Returns
A tuple `(fig, ax, states)` where:
- `fig`: The Makie `Figure`.
- `ax`: The axis containing the animated Wigner plot.
- `states`: The list of quantum states animated.
"""
function animate_wigner end




# ============================================================================ #

@testitem "Plot unitary populations" begin
    using CairoMakie
    using NamedTrajectories
    using PiccoloQuantumObjects

    H_drift = PAULIS[:X]
    H_drives = [PAULIS[:Z], PAULIS[:Y]]

    N = 100

    Δt = 0.1
    ts = collect(0:Δt:Δt*(N-1))

    a = 0.1 * randn(length(H_drives), length(ts))

    Us = exp.(-im * [(H_drift + sum(a[:, k] .* H_drives)) * ts[k] for k = 1:N])

    traj = NamedTrajectory(
        (
            Ũ⃗ = hcat(operator_to_iso_vec.(Us)...),
            a = a
        );
        controls = :a,
        timestep = Δt
    )

    fig = plot_unitary_populations(traj)
    save("../assets/unitary_populations.png", fig)
    @test fig isa Figure
end

end