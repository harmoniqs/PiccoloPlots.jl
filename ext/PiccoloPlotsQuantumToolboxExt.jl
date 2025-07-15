module PiccoloPlotsQuantumToolboxExt
using PiccoloPlots
using QuantumToolbox
using Makie
using NamedTrajectories
using PiccoloQuantumObjects
using TestItemRunner
using TestItems



function PiccoloPlots.plot_bloch(
    ::Val{:Makie},
    traj::NamedTrajectory{};
    state_name::Symbol = :ψ̃,
    kwargs...
)
    iso_states = eachcol(traj[state_name])
    kets = iso_to_ket.(iso_states)
    states = QuantumObject.(kets)
    bloch_vectors = QuantumToolbox._state_to_bloch.(states)

    b = QuantumToolbox.Bloch()
    QuantumToolbox.add_points!(b, hcat(bloch_vectors...))
    # QuantumToolbox.add_vectors!(b, hcat(bloch_vectors...)) # this doesn't work because it is a matrix and not a vector for hcat
    fig, lscene = QuantumToolbox.render(b; kwargs...)
    display(fig)
    return fig, lscene, states
end 

function PiccoloPlots.plot_bloch(
    traj::NamedTrajectory; kwargs...
)
    PiccoloPlots.plot_bloch(Val(:Makie), traj; kwargs...)
end

function PiccoloPlots.plot_wigner(
    ::Val{:Makie},
    traj::NamedTrajectory;
    state_name::Symbol = :ψ̃,
    xvec = -5:0.1:5,
    yvec = -5:0.1:5,
    projection::Val = Val(:two_dim),
    colorbar::Bool = true,
    kwargs...
)
    iso_states = eachcol(traj[state_name])
    kets = iso_to_ket.(iso_states)
    states = QuantumObject.(kets)

    fig, ax, hm = QuantumToolbox.plot_wigner(
        states[1];
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
    using PiccoloPlots
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
    println(traj)

    fig, lscene, states = PiccoloPlots.plot_bloch(Val(:Makie), traj)

    @test fig isa Figure
    @test lscene isa LScene
end

@testitem "plot_bloch shows expected curved Bloch path" begin
    using PiccoloPlots
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects: ket_to_iso
    using CairoMakie

    T = 20 # number of time steps
    ts = range(0, π/2; length=T)  # quarter-turn from |0⟩ to |1⟩

    #contstructing kets
    kets = [QuantumObject(cos(θ) * [1.0 + 0im, 0.0 + 0im] + sin(θ) * [0.0 + 0im, 1.0 + 0im]) for θ in ts]

    iso_kets = ket_to_iso.(ψ.data for ψ in kets)

    ψ̃ = hcat(iso_kets...)
    Δt = fill(1.0, T)

    comps = (
        ψ̃ = ψ̃,
        Δt = Δt,
    )

    traj = NamedTrajectory(comps)

    fig, lscene = PiccoloPlots.plot_bloch(Val(:Makie), traj)

    @test isa(fig, Figure)
    @test isa(lscene, LScene)
end

@testitem "plot_bloch with problem-constructed trajectory" begin
    using PiccoloPlots
    using QuantumToolbox
    using NamedTrajectories
    using PiccoloQuantumObjects
    using QuantumCollocation
    using CairoMakie

    CairoMakie.activate!()
    
    T = 50
    Δt = 0.2
    sys = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]])
    ψ_inits = Vector{ComplexF64}.([[1.0, 0.0], [0.0, 1.0]])
    ψ_targets = Vector{ComplexF64}.([[0.0, 1.0], [1.0, 0.0]])
    println("I am getting here")

    prob = QuantumStateSmoothPulseProblem(
        sys, ψ_inits, ψ_targets, T, Δt;
        piccolo_options=PiccoloOptions(verbose=false),
    )
    println("I am getting here")

    solve!(prob, max_iter=100, print_level=5)

    println("I am getting here")

    # traj = NamedTrajectory(comps)
    traj = prob.trajectory
    fig, lscene = PiccoloPlots.plot_bloch(Val(:Makie), traj, state_name=:ψ̃1)
    println("I am getting here")


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
        ψ̃ = hcat(ket_to_iso(ψ.data)),
        Δt = [1.0], 
    ))
    fig, ax, hm, _ = PiccoloPlots.plot_wigner(Val(:Makie), traj, state_name=:ψ̃)

    @test fig isa Figure
    @test ax isa Axis
    @test hm !== nothing
end

# @testitem "plot_wigner julia logo" begin
#     # This is a test that is copied from the QuantumToolbox.jl examples
#     using QuantumToolbox
#     using CairoMakie

#     CairoMakie.activate!(type = "svg", pt_per_unit = 1)
#     CairoMakie.enable_only_mime!(MIME"image/svg+xml"()) 
#     N = 30  # Cutoff of the Hilbert space for the harmonic oscillator
#     ρ = 2.5  # Amplitude of the coherent state
#     θ1 = π / 2
#     θ2 = π / 2 + 2π / 3
#     θ3 = π / 2 + 4π / 3
#     α1 = ρ * exp(im * θ1)
#     α2 = ρ * exp(im * θ2)
#     α3 = ρ * exp(im * θ3)
#     ψ = coherent(N, α1) + coherent(N, α2) + coherent(N, α3)
#     normalize!(ψ)
#     xvec = range(-ρ, ρ, 500) .* 1.5
#     yvec = xvec .+ (abs(imag(α1)) - abs(imag(α2))) / 2

#     fig = Figure(size = (250, 250), figure_padding = 0)
#     fig, ax, hm = plot_wigner(ψ, xvec = xvec, yvec = yvec, g = 2, library = Val(:Makie), location = fig[1,1])
#     hidespines!(ax)
#     hidexdecorations!(ax)
#     hideydecorations!(ax)
#     γ = 0.012

#     a = destroy(N)
#     H = a' * a
#     c_ops = [sqrt(γ) * a]

#     tlist = range(0, 2π, 100)

#     sol = mesolve(H, ψ, tlist, c_ops, progress_bar = Val(false))

#     fig = Figure(size = (250, 250), figure_padding = 0)
#     fig, ax, hm = plot_wigner(sol.states[end], xvec = xvec, yvec = yvec, g = 2, library = Val(:Makie), location = fig[1,1])
#     hidespines!(ax)
#     hidexdecorations!(ax)
#     hideydecorations!(ax)

#     fig, ax, hm

#     function set_color_julia(x, y, wig::T, α1, α2, α3, cmap1, cmap2, cmap3, δ) where {T}
#         amp1 = gaussian(x, real(α1), δ) * gaussian(y, imag(α1), δ)
#         amp2 = gaussian(x, real(α2), δ) * gaussian(y, imag(α2), δ)
#         amp3 = gaussian(x, real(α3), δ) * gaussian(y, imag(α3), δ)

#         c1 = get(cmap1, wig)
#         c2 = get(cmap2, wig)
#         c3 = get(cmap3, wig)

#         c_tot = (amp1 * c1 + amp2 * c2 + amp3 * c3) / (amp1 + amp2 + amp3)

#         wig_abs = abs(2 * (wig - 1 / 2))
#         # We introduce some non-linearity to increase the contrast
#         alpha = 2 * (1 / (1 + exp(-5 * wig_abs)) - 1 / 2)

#         return RGBAf(c_tot.r, c_tot.g, c_tot.b, alpha)
#     end

#     wig = wigner(sol.states[end], xvec, yvec, g = 2)
#     X, Y = meshgrid(xvec, yvec)
#     δ = 1.25 # Smoothing parameter for the Gaussian functions
#     julia_red = RGBAf(0.796, 0.235, 0.2, 1.0)
#     julia_green = RGBAf(0.22, 0.596, 0.149, 1.0)
#     julia_blue = RGBAf(0.251, 0.388, 0.847, 1.0)
#     julia_purple = RGBAf(0.584, 0.345, 0.698, 1.0)
#     n_repeats = 2

#     cmap1 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_green, n_repeats)))
#     cmap2 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_red, n_repeats)))
#     cmap3 = cgrad(vcat(fill(julia_blue, n_repeats), fill(julia_purple, n_repeats)))

#     vmax = maximum(wig)
#     wig_normalized = wig ./ (vmax * 2) .+ 1 / 2

#     img = set_color_julia.(X, Y, wig_normalized, α1, α2, α3, Ref(cmap1), Ref(cmap2), Ref(cmap3), δ)

#     fig = Figure(size = (250, 250), figure_padding = 0, backgroundcolor = :transparent)
#     ax = Axis(fig[1, 1], backgroundcolor = :transparent)
#     image!(ax, img', rasterize = 1)
#     hidespines!(ax)
#     hidexdecorations!(ax)
#     hideydecorations!(ax)
#     fig
# end


end