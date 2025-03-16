var documenterSearchIndex = {"docs":
[{"location":"lib/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"lib/","page":"Library","title":"Library","text":"Modules = [PiccoloPlots.QuantumObjectPlots]","category":"page"},{"location":"lib/#PiccoloPlots.QuantumObjectPlots.plot_unitary_populations-Tuple{NamedTrajectories.StructNamedTrajectory.NamedTrajectory}","page":"Library","title":"PiccoloPlots.QuantumObjectPlots.plot_unitary_populations","text":"plot_unitary_populations(\n    traj::NamedTrajectory;\n    unitary_columns::AbstractVector{Int}=1:2,\n    unitary_name::Symbol=:Ũ⃗,\n    control_name::Symbol=:a,\n    kwargs...\n)\n\nPlot the populations of the unitary columns of the unitary matrix in the trajectory. kwargs are passed to NamedTrajectories.plot.\n\nKeyword Arguments\n\nunitary_columns::AbstractVector{Int}: The columns of the unitary matrix to plot the populations of. Default is 1:2.\nunitary_name::Symbol: The name of the unitary matrix in the trajectory. Default is :Ũ⃗.\ncontrol_name::Symbol: The name of the control in the trajectory. Default is :a.\nkwargs...: Additional keyword arguments passed to NamedTrajectories.plot.\n\n\n\n\n\n","category":"method"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"EditURL = \"../../literate/quickstart.jl\"","category":"page"},{"location":"generated/quickstart/#Quickstart-Guide","page":"Quickstart Guide","title":"Quickstart Guide","text":"","category":"section"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Here is a simple example where we set up a NamedTrajectory with some dummy data and plot populations of the columns of the unitary matrix. First we will load some of the necessary packages:","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"using CairoMakie\nusing NamedTrajectories\nusing PiccoloQuantumObjects\nusing PiccoloPlots","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Next we will define some Hamiltonians","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"H_drives = [PAULIS[:Z], PAULIS[:Y]]\nH_drift = PAULIS[:X]","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Now we will generate some dummy control data","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"N = 100\nΔt = 0.1\nts = 0:Δt:Δt*(N-1)\nA = 0.1 * randn(length(H_drives), length(ts))","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Now we will generate the unitaries","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Us = exp.(-im * [(H_drift + sum(A[:, k] .* H_drives)) * ts[k] for k = 1:N])\nUs[1]","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"And create the trajectory","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"traj = NamedTrajectory(\n    (\n        Ũ⃗ = hcat(operator_to_iso_vec.(Us)...), # here we store the isomorphisms\n        a = A\n    );\n    controls = :a,\n    timestep = Δt\n)","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"Finally we will plot the populations","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"plot_unitary_populations(traj)","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"We can also only plot the first column (or any other subset of columns)","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"plot_unitary_populations(traj; unitary_columns=[1])","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"","category":"page"},{"location":"generated/quickstart/","page":"Quickstart Guide","title":"Quickstart Guide","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<div align=\"center\">\n  <a href=\"https://github.com/harmoniqs/Piccolo.jl\">\n    <img src=\"assets/logo.svg\" alt=\"logo\" width=\"25%\"/>\n  </a> \n</div>\n\n<div align=\"center\">\n  <table>\n    <tr>\n      <td align=\"center\">\n        <b>Documentation</b>\n        <br>\n        <a href=\"https://docs.harmoniqs.co/PiccoloPlots.jl/stable/\">\n          <img src=\"https://img.shields.io/badge/docs-stable-blue.svg\" alt=\"Stable\"/>\n        </a>\n        <a href=\"https://docs.harmoniqs.co/PiccoloPlots.jl/dev/\">\n          <img src=\"https://img.shields.io/badge/docs-dev-blue.svg\" alt=\"Dev\"/>\n        </a>\n      </td>\n      <td align=\"center\">\n        <b>Build Status</b>\n        <br>\n        <a href=\"https://github.com/harmoniqs/PiccoloPlots.jl/actions/workflows/CI.yml?query=branch%3Amain\">\n          <img src=\"https://github.com/harmoniqs/PiccoloPlots.jl/actions/workflows/CI.yml/badge.svg?branch=main\" alt=\"Build Status\"/>\n        </a>\n        <a href=\"https://codecov.io/gh/harmoniqs/PiccoloPlots.jl\">\n          <img src=\"https://codecov.io/gh/harmoniqs/PiccoloPlots.jl/branch/main/graph/badge.svg\" alt=\"Coverage\"/>\n        </a>\n      </td>\n      <td align=\"center\">\n        <b>License</b>\n        <br>\n        <a href=\"https://opensource.org/licenses/MIT\">\n          <img src=\"https://img.shields.io/badge/License-MIT-yellow.svg\" alt=\"MIT License\"/>\n        </a>\n      </td>\n      <td align=\"center\">\n        <b>Support</b>\n        <br>\n        <a href=\"https://unitary.fund\">\n          <img src=\"https://img.shields.io/badge/Supported%20By-Unitary%20Fund-FFFF00.svg\" alt=\"Unitary Fund\"/>\n        </a>\n      </td>\n    </tr>\n  </table>\n</div>\n\n<div align=\"center\">\n  <i>Easy plotting of quantum control trajectories</i>\n  <br>\n</div>","category":"page"},{"location":"#PiccoloPlots.jl","page":"Home","title":"PiccoloPlots.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PiccoloPlots.jl is designed to hold plotting recipes for the Piccolo.jl package. It utilzes the plot method of the NamedTrajectories.jl.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PiccoloPlots.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ] add PiccoloPlots","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As a simple example, we provide the plot_unitary_populations function which plots the populations of select columns of the unitary matrix. This can be employed as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using NamedTrajectories\nusing PiccoloQuantumObjects\nusing PiccoloPlots\n\n# Define the Hamiltonian H = X + a_1(t)Z + a_2(t)Y\nH_drift = PAULIS[:X]\nH_drives = [PAULIS[:Z], PAULIS[:Y]]\n\n# Generate control trajectory \nN = 100\n\nΔt = 0.1\nts = collect(0:Δt:Δt*(N-1))\n\na = 0.1 * randn(length(H_drives), length(ts))\n\n# Generate the unitaries\nUs = exp.(-im * [(H_drift + sum(a[:, k] .* H_drives)) * ts[k] for k = 1:N])\n\n# Create a NamedTrajectory\ntraj = NamedTrajectory(\n    (\n        Ũ⃗ = hcat(operator_to_iso_vec.(Us)...),\n        a = a\n    );\n    controls = :a,\n    timestep = Δt\n)\n\n# Plot the populations of the first and second qubits\nplot_unitary_populations(traj)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"}]
}
