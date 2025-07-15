using Test
using TestItemRunner
# using QuantumToolbox, Makie, GLMakie, Colors, PiccoloPlots
using QuantumCollocation, QuantumToolbox, PiccoloPlots # PiccoloPlotsQuantumToolboxExt must be loaded (via import of the associated weakdeps) to ensure the methods PiccoloPlots.plot_{bloch,wigner} are defined ahead of running package tests

# Run all testitem tests in package
@run_package_tests
