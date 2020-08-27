using Polymer
using Documenter

makedocs(;
    modules=[Polymer],
    authors="Yi-Xin Liu <lyx@fudan.edu.cn> and contributors",
    repo="https://github.com/liuyxpp/Polymer.jl/blob/{commit}{path}#L{line}",
    sitename="Polymer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://liuyxpp.github.io/Polymer.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/liuyxpp/Polymer.jl",
)
