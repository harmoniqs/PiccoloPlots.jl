module QuantumToolboxPlots

export animate_bloch
export animate_name
export animate_wigner
export plot_bloch!
export plot_wigner!

using TestItems

"""
    animate_bloch(traj::NamedTrajectory; state_name::Symbol=:ψ̃, fps::Int=60)

Animate the evolution of a quantum trajectory on the Bloch sphere without using Quantum Toolbox.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing quantum states.
- `state_name::Symbol`: The name of the quantum state in the trajectory. Default is `:ψ̃`.
- `fps::Int`: Frames per second for the animation. Default is `60`.

# Returns
The Makie `Figure` object displaying the animation.
"""
function animate_bloch end

function plot_bloch! end

"""
    animate_name(traj::NamedTrajectory; state_name::Symbol=:x, fps::Int=30)

Animate the evolution of a scalar or vector-valued variable in the trajectory.

# Arguments
- `traj::NamedTrajectory`: The trajectory containing the variable.
- `state_name::Symbol`: The name of the variable in the trajectory. Default is `:x`.
- `fps::Int`: Frames per second for the animation. Default is `30`.

# Returns
- `fig`: The Makie `Figure`.
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
- `fig`: The Makie `Figure`.
"""
function animate_wigner end

function plot_wigner! end


end