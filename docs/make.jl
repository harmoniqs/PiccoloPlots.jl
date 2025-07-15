using PiccoloPlots
using Documenter
using Literate

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

@info "Building Documenter site for PiccoloPlots.jl"
open(joinpath(@__DIR__, "src", "index.md"), write = true) do io
    for line in eachline(joinpath(@__DIR__, "..", "README.md"))
        if occursin("<!--", line) && occursin("-->", line)
            comment_content = match(r"<!--(.*)-->", line).captures[1]
            write(io, comment_content * "\n")
        else
            write(io, line * "\n")
        end
    end
end

pages = [
    "Home" => "index.md",
    "Quickstart Guide" => "generated/quickstart.md",
    "Library" => "lib.md"
]

format = Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://docs.harmoniqs.co/PiccoloPlots.jl", # TODO: should be "https://docs.harmoniqs.co/PiccoloPlots"; fix this in "github.com/harmoniqs/{DirectTrajOpt,NamedTrajectories,Piccolo,PiccoloPlots,PiccoloQuantumObjects,QuantumCollocation,TrajectoryIndexingUtils}.jl/tree/main/docs/make.jl"
    edit_link="main",
    assets=String[],
    mathengine = MathJax3(Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
            "tags" => "ams",
            "packages" => [
                "base",
                "ams",
                "autoload",
                "physics"
            ],
            "macros" => Dict(
                "minimize" => ["\\underset{#1}{\\operatorname{minimize}}", 1],
            )
        ),
        # :TeX => Dict(
        #     :Macros => Dict(
        #         :minimize => ["\\underset{#1}{\\operatorname{minimize}}", 1],
        #     )
        # )
    )),
)

src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "literate")
assets = joinpath(@__DIR__, "..", "assets")

lit_output = joinpath(src, "generated")
assets_output = joinpath(@__DIR__, "src", "assets")

for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>lit_output))[1]
    Literate.markdown(ipath, opath)
end

cp(assets, assets_output, force=true)

makedocs(;
    modules=[PiccoloPlots],
    authors="Aaron Trowbridge <aaron.j.trowbridge@gmail.com> and contributors",
    sitename="PiccoloPlots.jl",
    warnonly = [:missing_docs],
    format=format,
    pages=pages,
)

deploydocs(;
    repo="github.com/harmoniqs/PiccoloPlots.jl.git",
    devbranch="main",
)
