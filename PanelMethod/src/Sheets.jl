module Sheets

export Sheet, SheetGroup, aerodynamic_center

using Base.Iterators: drop
import LinearAlgebra
using LinearAlgebra: dot

using ..Utils
using ..Airfoils

"""
    Sheet{T<:AbstractFloat}

A connected sequence of line segments that define a vortex sheet.
"""
struct Sheet{T<:AbstractFloat}
    endpoints::Vector{Complex{T}}
    ctrl_points::Vector{Complex{T}}
    ctrl_normals::Vector{Complex{T}}
    @doc """
        Sheet{T}(npanels)
    
    Create an undefined sheet with `npanels` panels of float type `T`.
    """
    function Sheet{T}(npanels::Integer) where T
        npanels > 1 || throw(DomainError(npanels, "npanels must be greater than 1"))
    
        endpoints = Vector{Complex{T}}(undef, npanels + 1)
        ctrl_points = Vector{Complex{T}}(undef, npanels)
        ctrl_normals = Vector{Complex{T}}(undef, npanels)
        return new(endpoints, ctrl_points, ctrl_normals)
    end
end

"""
    Sheet(airfoil, npanels, offset=1e-8)

Create a sheet with `npanels` panels and vertices on `airfoil`.

`offset` determines the distance (as a fraction of chord length) that the control points
used for calculations are offset from each panel's midpoint.
"""
function Sheet(airfoil::Airfoil{T}, npanels::Integer) where T
    update!(Sheet{T}(npanels), airfoil)
end

function Sheet(airfoil::Airfoil{T}, npanels::Integer, offset::Real) where T
    update!(Sheet{T}(npanels), airfoil, offset)
end

"""
    update!(sheet::Sheet, airfoil::Airfoil, offset=1e-8)

The in place version of `Sheet(airfoil, n, offset)` if `sheet` has `n` panels.
"""
function update!(sheet::Sheet{T}, airfoil::Airfoil{T}, offset::Real = T(1e-8)) where T
    endpoints = sheet.endpoints
    s = LinRange{T}(0, 1, length(sheet.endpoints))
    endpoints .= Airfoils.position.(Ref(airfoil), s)
    
    rot = ifelse(is_clockwise(airfoil), 1im, -1im)
    offset_dist = chordlength(airfoil) * offset

    for (i, z1, z2) in zip(eachindex(sheet.ctrl_points), endpoints, drop(endpoints, 1))
        normal = normalize(rot * (z2 - z1))
        sheet.ctrl_points[i] = (z1 + z2) / 2 + offset_dist * normal
        sheet.ctrl_normals[i] = normal
    end

    return sheet
end

function _lengths(sheet::Sheet)
    # Return a generator over the lengths of each panel in `sheet`
    (abs(z1 - z2) for (z1, z2) in zip(sheet.endpoints, drop(sheet.endpoints, 1)))
end


"""
    SheetGroup{T<:AbstractFloat, S<:AbstractVector{Sheet{T}}}

A group of `Sheet`s and the solution to the linear system determining the vortex strength
distribution.

If any of the underlying sheets of `group` are mutated, `update!(group)` must be called
to resolve the linear system.
"""
struct SheetGroup{T<:AbstractFloat, S<:AbstractVector{Sheet{T}}}
    sheets::S
    _A::LinearAlgebra.Transpose{T, Matrix{T}}
    _x_re::Vector{T}
    _x_im::Vector{T}
    @doc """
        SheetGroup(sheets, [undef]))
    
    Create a `SheetGroup` from a single or vector of `Sheet`. If `undef` is provided,
    do not initialize the data with `update!`.
    """
    function SheetGroup(sheets::S, ::UndefInitializer
                        ) where {T, S<:AbstractVector{Sheet{T}}}
        γ_count = sum(_point_count, sheets)
        
        # since A is assigned to by rows, work on a transposed (or row-major) matrix
        A = transpose(Matrix{T}(undef, γ_count, γ_count))
        x_re = Vector{T}(undef, γ_count)
        x_im = Vector{T}(undef, γ_count)
        
        return new{T, S}(sheets, A, x_re, x_im)
    end
end

SheetGroup(sheets::AbstractVector{Sheet{T}}) where {T} = update!(SheetGroup(sheets, undef))

function SheetGroup(sheet::Sheet, args::Vararg{UndefInitializer, N}) where N
    SheetGroup([sheet], args...)
end

"""
    update!(group::SheetGroup)

Updates the data in `group`. Should be called if any of the underlying `Sheet`s are mutated.
"""
function update!(group::SheetGroup)
    A = _fill_A_matrix!(group._A, group.sheets)
    x_reim = _fill_b_vector!(group._x_re, group._x_im, group.sheets)

    qrA = LinearAlgebra.qr!(A)
    foreach(x -> LinearAlgebra.ldiv!(qrA, x), x_reim)
    return group
end

Base.iterate(group::SheetGroup) = iterate(group.sheets)
Base.iterate(group::SheetGroup, state) = iterate(group.sheets, state)
Base.length(group::SheetGroup) = length(group.sheets)
Base.eltype(::Type{SheetGroup{T}}) where {T} = Sheet{T}
Base.IndexStyle(::Type{SheetGroup}) = IndexLinear()


_point_count(sheet::Sheet) = length(sheet.endpoints)
_point_count(group::SheetGroup) = sum(_point_count, group)


"""
    aerodynamic_center(sheets, chord::AbstractChord)

Returns the aerodynamic center `t` of `sheets` as a fraction of chord length.

`chordline(chord, t)` retrieves the coordinate of the aerodynamic center. `sheets`
may be of type `Sheet`, `AbstractVector{Sheet}`, or `SheetGroup`.
"""
function aerodynamic_center(sheetgroup::SheetGroup{T}, chord::AbstractChord) where T
    c::Complex{T} = chordvec(chord)
    conj_unit_c::Complex{T} = conj(normalize(c))
    LE::Complex{T} = leading_edge(chord)
    
    num::T = 0
    den::T = 0
    x_re, x_im = sheetgroup._x_re, sheetgroup._x_im
    for sheet in sheetgroup
        for (L, pt, normal) in zip(_lengths(sheet), sheet.ctrl_points, sheet.ctrl_normals)
            tangent = im * normal

            indexes = Iterators.Stateful(eachindex(x_re))
            d::Complex{T} = 0
            for src_sheet in sheetgroup
                local a::T
                next_a::T = 0
                for (z1, z2) in zip(src_sheet.endpoints, drop(src_sheet.endpoints, 1))
                    a1, a2 = cdot.(tangent, _γ_coefficients(pt, z1, z2))
                    a = next_a + a1
                    next_a = a2

                    i = popfirst!(indexes)
                    d += a * complex(x_re[i], x_im[i])
                end
                i = popfirst!(indexes)
                d += next_a * complex(x_re[i], x_im[i])
            end

            tmp1 = conj_unit_c * (tangent + d)
            tmp2 = L * real(tmp1) * imag(tmp1)
            num += tmp2 * ccross(pt - LE, normal)
            den += tmp2 * ccross(c, normal)
        end
    end

    return num / den
end

function aerodynamic_center(sheets, chord::AbstractChord) where T
    aerodynamic_center(SheetGroup(sheets), chord)
end

function _fill_A_matrix!(A::AbstractMatrix{T}, sheets::AbstractVector{Sheet{T}}) where T
    next_i = 1
    for sheet in sheets
        i = next_i
        next_i += _point_count(sheet)

        _γ_coefficients!(
            view(A, i:next_i-2, :), sheet.ctrl_points, sheet.ctrl_normals, sheets)

        kutta_i = next_i - 1
        fill!(view(A, kutta_i, :), 0)
        A[kutta_i, i] = 1
        A[kutta_i, kutta_i] = 1
    end
    return A
end

function _γ_coefficients!(out::AbstractMatrix{T},
                          points::Vector{Complex{T}},
                          alongs::Vector{Complex{T}},
                          sheets::AbstractVector{Sheet{T}}) where T
    # Fill each row in `out` with the vortex/source strength coefficients of `sheets` for
    # each point and along vector in `points` and `alongs`.
    for (row, point, along) in zip(eachrow(out), points, alongs)
        _γ_coefficients!(row, point, along, sheets)
    end
end

function _γ_coefficients!(out::AbstractVector{T},
                          point::Complex{T},
                          along::Complex{T},
                          sheets::AbstractVector{Sheet{T}}) where T
    indices = eachindex(out)
    i = 1
    for sheet in sheets
        this_i = i
        i += _point_count(sheet)
        _γ_coefficients!(view(out, indices[this_i:i-1]), point, along, sheet)
    end
end

function _γ_coefficients!(out::AbstractVector{T},
                          point::Complex{T},
                          along::Complex{T},
                          sheet::Sheet{T}) where T
    indexes = eachindex(out)
    out[begin] = zero(T)
    @inbounds for (i1, i2, z1, z2) in
            zip(indexes, drop(indexes, 1), sheet.endpoints, drop(sheet.endpoints, 1))
        a1, a2 = _γ_coefficients(point, z1, z2)
        out[i1] += cdot(a1, along)
        out[i2] = cdot(a2, along)
    end
end

function _γ_coefficients(z::Complex{T}, z1::Complex{T}, z2::Complex{T}) where T
    # z1 and z2 define the endpoints of a linearly distributed vortex panel with
    # endpoint vortex strengths γ1 and γ2. Return the coefficients of γ1 and γ2
    # in the flow velocity at z.
    tmp1 = z1 - z2
    tmp2 = -abs(tmp1) / (2im * pi * tmp1^2)
    tmp3 = tmp2 * log((z - z1) / (z - z2))
    tmp4 = tmp1 * tmp2
    return conj.((tmp3 * (z2 - z) - tmp4, tmp3 * (z - z1) + tmp4))
end

function _fill_b_vector!(b_re::AbstractVector{T}, b_im::AbstractVector{T},
                         sheets::AbstractVector{Sheet{T}}) where T
    @assert axes(b_re) == axes(b_im)

    indexes = Iterators.Stateful(eachindex(b_re))
    for sheet in sheets
        for (i, normal) in zip(indexes, sheet.ctrl_normals)
            b_re[i], b_im[i] = reim(-normal)
        end
        i = popfirst!(indexes)
        b_re[i] = b_im[i] = 0
    end
    return (b_re, b_im)
end

end # module
