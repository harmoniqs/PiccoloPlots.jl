using PiccoloPlots
using PiccoloDocsTemplate

pages = [
    "Home" => "index.md",
    "Quickstart Guide" => "generated/quickstart.md",
    "Library" => "lib.md"
]

generate_docs(
    @__DIR__,
    "PiccoloPlots",
    PiccoloPlots,
    pages;
    format_kwargs = (canonical = "https://docs.harmoniqs.co/PiccoloPlots.jl",),
)