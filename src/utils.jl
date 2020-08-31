"""
    unicodesymbol2string(us::Char)

Return an ascii representation (LaTeX name for a unicode symbol) for the input character. If the input character is already an ASCII or it is not listed in `REPL.REPLCompletions.latex_symbols`, the input string is returned.

```julia-repl
julia> unicodesymbol2string("β")
"beta"
julia> unicodesymbol2string("b")
"b"
```
"""
function unicodesymbol2string(us::Char)
    us = string(us)
    return isascii(us) ? us : isempty(symbol_latex(us)) ? us : symbol_latex(us)[2:end]
end