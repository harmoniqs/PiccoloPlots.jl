<!--```@raw html-->
<div align="center">
  <a href="https://github.com/harmoniqs/Piccolo.jl">
    <img src="assets/logo.svg" alt="logo" width="25%"/>
  </a> 
</div>

<div align="center">
  <table>
    <tr>
      <td align="center">
        <b>Documentation</b>
        <br>
        <a href="https://docs.harmoniqs.co/PiccoloPlots/stable/">
          <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"/>
        </a>
        <a href="https://docs.harmoniqs.co/PiccoloPlots/dev/">
          <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"/>
        </a>
      </td>
      <td align="center">
        <b>Build Status</b>
        <br>
        <a href="https://github.com/harmoniqs/PiccoloPlots.jl/actions/workflows/CI.yml?query=branch%3Amain">
          <img src="https://github.com/harmoniqs/PiccoloPlots.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"/>
        </a>
        <a href="https://codecov.io/gh/harmoniqs/PiccoloPlots.jl">
          <img src="https://codecov.io/gh/harmoniqs/PiccoloPlots.jl/branch/main/graph/badge.svg" alt="Coverage"/>
        </a>
      </td>
      <td align="center">
        <b>License</b>
        <br>
        <a href="https://opensource.org/licenses/MIT">
          <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
        </a>
      </td>
      <td align="center">
        <b>Support</b>
        <br>
        <a href="https://unitary.fund">
          <img src="https://img.shields.io/badge/Supported%20By-Unitary%20Fund-FFFF00.svg" alt="Unitary Fund"/>
        </a>
      </td>
    </tr>
  </table>
</div>

<div align="center">
  <i>Easy plotting of quantum control trajectories</i>
  <br>
</div>
<!--```-->

# PiccoloPlots.jl

**PiccoloPlots.jl** is designed to hold plotting recipes for the Piccolo.jl package. It utilzes the `plot` method of the [`NamedTrajectories.jl`](https://github.com/harmoniqs/NamedTrajectories.jl).

### Installation
PiccoloPlots.jl can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
julia> ] add PiccoloPlots
```

### Usage
As a simple example, we provide the `plot_unitary_populations` function which plots the populations of select columns of the unitary matrix. This can be employed as follows:

```julia
using NamedTrajectories
using PiccoloQuantumObjects
using PiccoloPlots

# Define the Hamiltonian H = X + u_1(t)Z + u_2(t)Y
H_drift = PAULIS[:X]
H_drives = [PAULIS[:Z], PAULIS[:Y]]

# Generate control trajectory 
N = 100

Δt = 0.1
ts = collect(0:Δt:Δt*(N-1))

u = 0.1 * randn(length(H_drives), length(ts))

# Generate the unitaries
Us = exp.(-im * [(H_drift + sum(u[:, k] .* H_drives)) * ts[k] for k = 1:N])

# Create a NamedTrajectory
traj = NamedTrajectory(
    (
        Ũ⃗ = hcat(operator_to_iso_vec.(Us)...),
        u = u,
        Δt = fill(Δt, N),
    );
    controls = :u,
    timestep = :Δt,
)

# Plot the populations of the first and second qubits
plot_unitary_populations(traj)
```
```

![](assets/unitary_populations.png)


### Building Documentation
This package uses a Documenter config that is shared with many of our other repositories. To build the docs, you will need to run the docs setup script to clone and pull down the utility. 
```
# first time only
./docs/get_docs_utils.sh   # or ./get_docs_utils.sh if cwd is in ./docs/
```

To build the docs pages:
```
julia --project=docs docs/make.jl
```

or editing the docs live:
```
julia --project=docs
> using LiveServer, PiccoloPlots, Revise
> servedocs(literate_dir="docs/literate", skip_dirs=["docs/src/generated", "docs/src/assets/"], skip_files=["docs/src/index.md"])
```

> **Note:** `servedocs` needs to watch a subset of the files in the `docs/` folder. If it watches files that are generated on a docs build/re-build, `servedocs` will continuously try to re-serve the pages.
> 
> To prevent this, ensure all generated files are included in the skip dirs or skip files args for `servedocs`.

For example, if we forget index.md like so:
```
julia --project=docs
> using LiveServer, PiccoloPlots, Revise
> servedocs(literate_dir="docs/literate", skip_dirs=["docs/src/generated", "docs/src/assets/"])
```
it will not build and serve.

-----

*"It seems that perfection is attained not when there is nothing more to add, but when there is nothing more to take away." - Antoine de Saint-Exupéry*