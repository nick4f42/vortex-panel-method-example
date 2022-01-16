module AirfoilTools

export airfoil_names, get_airfoil

import Downloads
using InlineStrings: String31

using ..Airfoils: InterpolatedAirfoil


const SITE_BASE = "http://airfoiltools.com"
const SEARCH_URL = SITE_BASE * "/search/airfoils"
const DATFILE_URL = SITE_BASE * "/airfoil/lednicerdatfile"

const StringN = String31  # expect that airfoil names from airfoiltools.com are < 32 bytes
const _name_cache = StringN[]  # cache result from airfoil_names()
const NAMES_COUNT_HINT = 1750  # expected length of _name_cache

function airfoil_names()
    length(_name_cache) > 0 && return _name_cache
    
    sizehint!(_name_cache, NAMES_COUNT_HINT)

    html_str = (String ∘ take!)(Downloads.download(SEARCH_URL, IOBuffer()))
    for href in _each_href(html_str)
        m = match(r"airfoil=(.+)$", href)
        isnothing(m) && continue
        s = m[1]
        if sizeof(s) < sizeof(StringN)
            push!(_name_cache, s)
        else
            @warn "$(repr(s)) is ≥ $(sizeof(StringN)) bytes and will be ignored."
        end
    end
    return _name_cache
end

# captures the href property of an <a> tag in the "href" group
const _HREF_REGEX = r"<a[^>]*href\s*=\s*([\"'])(?<href>(?:(?!(?<!\\)\1).)*(?:(?!(?<!\\)\1).)*)\1"

function _each_href(html_str::AbstractString)
    # Return a generator over each href of each <a href="..."> tag.
    return (replace(m[:href], "\\"*m[1] => m[1])
            for m in eachmatch(_HREF_REGEX, html_str))
end


function get_airfoil(name::AbstractString)
    points, LE_idx = get_airfoil_points(name)
    return InterpolatedAirfoil(points, points[LE_idx], clockwise=false)
end

function get_airfoil_points(name::AbstractString)
    url = DATFILE_URL * "?airfoil=" * name
    str = (String ∘ take!)(Downloads.download(url, IOBuffer()))

    coord_iter = (
        map(x -> parse(Float64, x), split(line))
        for line in split(str, '\n', keepempty=false)
        if !all(isspace, line) && !occursin(r"[^\-.\d\s]", line)
    )
    
    first_nums, rest = Iterators.peel(coord_iter)
    if all(>(1), first_nums)
        upper, lower = map(Int, first_nums)
        coords = map(Base.splat(complex), rest)
        return _get_lednicer_points(coords, upper, lower)
    else
        coords = map(Base.splat(complex), coord_iter)
        return _get_selig_points(coords)
    end
end

function _get_lednicer_points(coords::Vector{T}, upper::Int, lower::Int) where T <: Complex
    @assert length(coords) == upper + lower

    points = Vector{T}(undef, length(coords) - 1)
    points[1:upper] .= coords[upper:-1:begin]
    points[upper+1:end] .= coords[upper+2:end]

    return (points, upper)
end

function _get_selig_points(coords::Vector{Complex})
    LE_idx = argmin(abs2 ∘ real, coords)
    return (coords, LE_idx)
end

end # module
