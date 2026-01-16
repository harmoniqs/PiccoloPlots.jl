module QuantumObjectPlots

export plot_unitary_populations
export plot_state_populations

using LaTeXStrings
using Makie
import Makie: plot
using LinearAlgebra
using NamedTrajectories
using PiccoloQuantumObjects
using TestItems

# Utility functions
function get_layout(index::Int, n::Int)
    √n = isqrt(n) + 1
    return ((index - 1) ÷ √n + 1, (index - 1) % √n + 1)
end

# Include submodules
include("unitary_populations.jl")
include("state_populations.jl")

end
