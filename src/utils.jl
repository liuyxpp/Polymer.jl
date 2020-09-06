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

"""
    reverse_dict(d::Dict)

Return a Dict instance with each `key => val` pair in `d` reversed as `val => key` pair. Same values in `d` will be overwirten by the last occurd pair.
"""
reverse_dict(d::Dict) = Dict(v => k for (k, v) in d)

"""
    _sort_tuple2(t)

Sort two-element Tuple to increase order. Example

```julia-REPL
julia> _sort_tuple2((2, 1))
(1, 2)
julia> _sort_tuple2((1, 2))
(1, 2)
```
"""
function _sort_tuple2(t)
	a, b = t
	return a > b ? (b, a) : (a, b)
end