# PiccoloPlots.jl Context

> AI-friendly context for maintaining consistency. Update this when making significant changes.

## Package Purpose

PiccoloPlots.jl provides **visualization utilities** for the Piccolo quantum control ecosystem. It builds on NamedTrajectories plotting and adds quantum-specific visualizations.

## Core Features

### Trajectory Plotting

Plot quantum control trajectories with component selection:

```julia
using PiccoloPlots
using NamedTrajectories

# Plot specific components
plot(traj, [:u, :Ũ⃗])

# Plot all components
plot(traj)

# Save to file
plot("output.png", traj, [:u])

# With options
plot(traj, [:u]; 
    title="Optimal Controls",
    xlabel="Time (ns)",
    ylabel="Amplitude"
)
```

### Animated Plots

Create animations of optimization progress:

```julia
using PiccoloPlots.AnimatedPlots

# During optimization, collect trajectory history
callback, trajectory_history = trajectory_history_callback(prob)
solve!(prob; callback=callback)

# Create animation frames
for (iter, traj) in enumerate(trajectory_history)
    plot("frame_$iter.png", traj, [:u])
end
```

### Quantum Object Visualization

Visualize quantum states and operators:

```julia
using PiccoloPlots.QuantumObjectPlots

# Plot Wigner function
plot_wigner(ρ)

# Plot state populations
plot_populations(ψ)

# Plot Bloch sphere
plot_bloch(ψ)
```

### QuantumToolbox Integration

Plots for QuantumToolbox.jl objects:

```julia
using PiccoloPlots.QuantumToolboxPlots
using QuantumToolbox

# Plot quantum state from QuantumToolbox
ket = basis(2, 0)
plot_state(ket)
```

## Module Structure

```
PiccoloPlots.jl/src/
├── PiccoloPlots.jl           # Main module
├── animated_plots.jl         # AnimatedPlots module
├── quantum_object_plots/
│   └── _quantum_object_plots.jl  # QuantumObjectPlots module
└── quantum_toolbox_plots.jl  # QuantumToolboxPlots module
```

## Integration with Other Packages

### NamedTrajectories

PiccoloPlots extends the base `plot` function from NamedTrajectories:

```julia
# Base NamedTrajectories plot
using NamedTrajectories
plot(traj)  # Basic plot

# Enhanced with PiccoloPlots
using PiccoloPlots
plot(traj)  # Quantum-aware styling
```

### QuantumCollocation

Plot results from quantum control problems:

```julia
using QuantumCollocation
using PiccoloPlots

# Solve problem
qcp = SmoothPulseProblem(qtraj, N)
solve!(qcp)

# Plot results
traj = get_trajectory(qcp)
plot(traj, [:u])           # Controls
plot(traj, [:Ũ⃗])          # Unitary evolution
plot(traj, [:u, :du])      # Controls and derivatives
```

## Common Plot Types

| Component | Description | Use |
|-----------|-------------|-----|
| `:u` | Control pulses | Visualize optimized drives |
| `:Ũ⃗` | Isomorphic unitary | Track unitary evolution |
| `:ψ̃` | Isomorphic ket | Track state evolution |
| `:du`, `:ddu` | Derivatives | Check smoothness |
| `:Δt` | Timesteps | Variable time discretization |

## Customization

### Plot Options

```julia
plot(traj, [:u];
    # Figure
    size=(800, 600),
    dpi=150,
    
    # Axes
    xlims=(0, T),
    ylims=(-1, 1),
    xlabel="Time",
    ylabel="Amplitude",
    
    # Style
    linewidth=2,
    legend=:topright,
    
    # Layout
    layout=(2, 1)  # Stack components
)
```

### Component-Specific Limits

```julia
plot(traj, [:u, :Ũ⃗];
    ylims=(u = (-1, 1), Ũ⃗ = (-2, 2))
)
```

## Dependencies

- **Plots.jl** - Base plotting backend
- **NamedTrajectories.jl** - Trajectory data structure
- **QuantumToolbox.jl** (optional) - Quantum state objects

## Typical Workflow

```julia
using Piccolo  # Includes PiccoloPlots

# Define and solve problem
sys = TransmonSystem([1.0, 1.0])
qtraj = UnitaryTrajectory(sys, GATES[:H], 10.0)
qcp = SmoothPulseProblem(qtraj, 51)
solve!(qcp)

# Visualize
traj = get_trajectory(qcp)

# Quick view
plot(traj)

# Publication-quality
plot("figure1.pdf", traj, [:u]; 
    size=(400, 300),
    xlabel="Time (ns)",
    ylabel="Drive (GHz)"
)
```

## Animation Example

```julia
using PiccoloPlots

# Collect history during optimization
callback, history = trajectory_history_callback(prob)
solve!(prob; callback=callback, max_iter=100)

# Create GIF frames
for (i, traj) in enumerate(history)
    idx = lpad(i, 3, '0')
    plot("frame_$idx.png", traj, [:u]; 
        ylims=(-1, 1),
        title="Iteration $i"
    )
end

# Combine with external tool (e.g., ffmpeg, ImageMagick)
# convert -delay 10 frame_*.png animation.gif
```
