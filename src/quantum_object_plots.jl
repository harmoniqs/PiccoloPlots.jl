module QuantumObjectPlots

export plot_unitary_populations
export plot_bloch

using LaTeXStrings
# only need MakieCore for plot name
using MakieCore
using LinearAlgebra
using GLMakie
using GeometryBasics
using NamedTrajectories
using PiccoloQuantumObjects
using QuantumCollocation
using TestItemRunner

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
    plot_bloch_trajectory(
        traj::NamedTrajectory;
        state_name::Symbol = :ψ̃,
        kwargs...
    )

Plot the Bloch sphere trajectory of a qubit state over time, with an interactive slider showing progression.

# Keyword Arguments
- `state_name::Symbol`: The name of the state vector in the trajectory. Default is `:ψ̃`.
- `kwargs...`: Additional keyword arguments passed to [`NamedTrajectories.plot`](https://docs.harmoniqs.co/NamedTrajectories/dev/generated/plotting/).
"""
function plot_bloch(
    traj::NamedTrajectory;
    state_name::Symbol = :ψ̃,
    kwargs...
)
    iso_vecs = eachcol(traj[state_name])
    kets = iso_to_ket.(iso_vecs)
    bloch_vectors = ket_to_bloch.(kets)

    bx = getindex.(bloch_vectors, 1)
    by = getindex.(bloch_vectors, 2)
    bz = getindex.(bloch_vectors, 3)
    points = Point3f.(bx, by, bz)

    θ = LinRange(0, π, 50)
    ϕ = LinRange(0, 2π, 50)
    x = [sin(t) * cos(p) for t in θ, p in ϕ]
    y = [sin(t) * sin(p) for t in θ, p in ϕ]
    z = [cos(t) for t in θ, p in ϕ]

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:equal, title="Bloch Sphere")
    wireframe!(ax, x, y, z, color=:lightgray, transparency=true)
    lines!(ax, points; color=:blue, linewidth=3, kwargs...)

    index = Observable(1)
    current_point = @lift Point3f(bx[$index], by[$index], bz[$index])
    scatter!(ax, current_point, color=:red, markersize=10)

    slider = Slider(fig[2, 1], range=1:length(points), startvalue=1)
    connect!(index, slider.value)

    return fig
end


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

@testitem "Plot on bloch" begin
    using CairoMakie
    using NamedTrajectories
    using PiccoloQuantumObjects

    θs = LinRange(0, 2π, 100)
    ψs = [cos(θ/2)*ket(0) + sin(θ/2)*im*ket(1) for θ in θs]
    ψs_mat = hcat(ket_to_iso.(ψs)...)

    traj = NamedTrajectory((ψ̃ = ψs_mat, ); timestep=0.1)

    fig = plot_blochy(traj)
    save("../assets/bloch_plot_w_slider.png", fig)
    @test fig isa Figure
end
end
