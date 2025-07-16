module QuantumObjectPlots

export plot_unitary_populations

using LaTeXStrings
# only need MakieCore for plot
using MakieCore
using LinearAlgebra
# using GeometryBasics
using NamedTrajectories
using PiccoloQuantumObjects
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
            a = a,
            Δt = ts,
        );
        controls = :a,
        timestep = :Δt,
    )

    fig = plot_unitary_populations(traj)
    save("../assets/unitary_populations.png", fig)
    @test fig isa Figure
end

end