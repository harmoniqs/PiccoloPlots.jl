module PiccoloPlots

using Reexport

# Ensure plotting dependencies are loaded
@reexport using NamedTrajectories: plot
using CairoMakie

include("animated_plots.jl")
@reexport using .AnimatedPlots

include("quantum_object_plots/_quantum_object_plots.jl")
@reexport using .QuantumObjectPlots

include("quantum_toolbox_plots.jl")
@reexport using .QuantumToolboxPlots

end
