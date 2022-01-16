module Airfoils

export
    AbstractChord, leading_edge, trailing_edge, chordvec, chordunit, chordlength, chordline,
    Airfoil, is_clockwise, Xform, mirror, JoukowskyAirfoil, InterpolatedAirfoil

using Base.Iterators: drop
using Interpolations

using ..Utils


"""
    AbstractChord

Supertype for types that have some sense of a "chord" with leading and trailing edges.
"""
abstract type AbstractChord end

"""
    leading_edge(chord)

The complex position of the `chord`'s leading edge.
"""
function leading_edge end

"""
    trailing_edge(chord)

The complex position of the `chord`'s trailing edge.
"""
function trailing_edge end

"""
    chordvec(chord)

Shorthand for `trailing_edge(chord) - leading_edge(chord)`.
"""
chordvec(chord::AbstractChord) = trailing_edge(chord) - leading_edge(chord)

"""
    chordunit(chord)

Unit complex number in the direction of the chord line. The normalized form of [`chordvec`](@ref).
"""
chordunit(chord::AbstractChord) = normalize(chordvec(chord))

"""
    chordlength(chord)

The chord length (absolute value of [`chordvec`](@ref)).
"""
chordlength(chord::AbstractChord) = abs(chordvec(chord))

"""
    chordline(chord, t)

Return a complex point a fraction of `t` from the leading to trailing edge.
"""
function chordline(chord::AbstractChord, t::Real)
    LE = leading_edge(chord)
    TE = trailing_edge(chord)
    return LE + t * (TE - LE)
end


"""
    Airfoil{T<:AbstractFloat} <: AbstractChord

Supertype for airfoil parameterizations.
"""
abstract type Airfoil{T<:AbstractFloat} <: AbstractChord end

"""
    position(airfoil, t)

Complex coordinate of `airfoil`'s parameterization at `t`, where `0 ≤ t ≤ 1`.

`position(airfoil, 0) == position(airfoil, 1) == trailing_edge(airfoil)`
"""
function position end

"""
    is_clockwise(airfoil)

Returns a `Bool` indicating whether `airfoil`'s parameterization is clockwise.
"""
function is_clockwise end


"""
    Xform{T} <: AbstractChord

A chord specified by the midpoint and difference of the leading and trailing edge.
"""
struct Xform{T} <: AbstractChord
    center::Complex{T}
    scale::Complex{T}
end

"""
    Xform(center, scale)

A chord specified by the midpoint `center` and chord vector `scale`.
"""
Xform(center::Number, scale::Number) = Xform(complex.(promote(center, scale))...)

"""
    Xform(airfoil)

Retrieve the equivalent [`Xform`](@ref) of `airfoil`.
"""
function Xform(airfoil::Airfoil{T}) where T
    LE = leading_edge(airfoil)
    TE = trailing_edge(airfoil)
    return Xform((LE + TE) / 2, TE - LE)
end

leading_edge(xform::Xform) = xform.center - xform.scale / 2
trailing_edge(xform::Xform) = xform.center + xform.scale / 2

Base.convert(t::Type{<:Xform}, xform::Xform) = t(xform.center, xform.scale)


"""
    XformAirfoil{T, A<:Airfoil{T}} <: Airfoil{T}

A transformed version of a base airfoil.
"""
struct XformAirfoil{T, A<:Airfoil{T}} <: Airfoil{T}
    base::A
    xform::Xform{T}
    base_xform::Xform{T}
    @doc """
        XformAirfoil(base, xform)

    A transformed version of `base` to align with the chord of `xform`.
    """
    function XformAirfoil(base::A, xform::Xform) where {T, A<:Airfoil{T}}
        return new{T, A}(base, xform, Xform(base))
    end
end

XformAirfoil(airfoil::XformAirfoil, xform::Xform) = XformAirfoil(airfoil.base, xform)

Xform(airfoil::XformAirfoil) = airfoil.xform

"""
    (xform::Xform)(airfoil)

Transforms `airfoil` to the chord of `xform`. Shorthand for `XformAirfoil(airfoil, xform)`. 
"""
(xform::Xform)(airfoil::Airfoil) = XformAirfoil(airfoil, xform)

leading_edge(airfoil::XformAirfoil) = leading_edge(airfoil.xform)
trailing_edge(airfoil::XformAirfoil) = trailing_edge(airfoil.xform)
is_clockwise(airfoil::XformAirfoil) = is_clockwise(airfoil.base)

function position(airfoil::XformAirfoil, t::Real)
    xf = airfoil.xform
    base_xf = airfoil.base_xform
    z = position(airfoil.base, t)
    return xf.center + xf.scale / base_xf.scale * (z - base_xf.center)
end


"""
    MirrorAirfoil{T, A<:Airfoil{T}} <: Airfoil{T}

A mirrored version of an airfoil.
"""
struct MirrorAirfoil{T, A<:Airfoil{T}} <: Airfoil{T}
    base::A
end

"""
    MirrorAirfoil(base)

A mirrored version of `base` across its chord line.
"""
MirrorAirfoil(airfoil::MirrorAirfoil) = airfoil.base

"""
    mirror(base)

A mirrored version of `base` across its chord line. Alias for [`MirrorAirfoil`](@ref).
"""
const mirror = MirrorAirfoil

for func in (:leading_edge, :trailing_edge, :Xform)
    @eval $func(airfoil::MirrorAirfoil) = $func(airfoil.base)
end

is_clockwise(airfoil::MirrorAirfoil) = !is_clockwise(airfoil.base)

function position(airfoil::MirrorAirfoil, t::Real)
    xform = Xform(airfoil.base)
    z = position(airfoil.base, t)
    return xform.center + xform.scale * conj((z - xform.center) / xform.scale)
end


"""
    JoukowskyAirfoil{T} <: Airfoil{T}

An airfoil defined using the Joukowsky transformation.
"""
struct JoukowskyAirfoil{T} <: Airfoil{T}
    μ::Complex{T}
    _LE_param::T  # parameter of the leading edge
    @doc """
        JoukowskyAirfoil(thickness, camber)

    `thickness` and `camber` define the real and imaginary parts of the circle's center
    before the Joukowsky transformation.
    """
    function JoukowskyAirfoil(thickness::T, camber::T) where T<:AbstractFloat
        if !(0 <= thickness <= 0.5 && 0 <= camber <= 0.5)
            throw(DomainError("thickness and camber must be non-negative."))
        end
        μ = complex(-thickness, camber)
        LE_param = mod(atan((real(μ) - 1) / imag(μ)) / pi, 1)
        return new{T}(μ, LE_param)
    end
end

function JoukowskyAirfoil{T}(thickness::Real, camber::Real) where T
    JoukowskyAirfoil(convert(T, thickness), convert(T, camber))
end

function JoukowskyAirfoil(thickness::T, camber::T) where T<:Real
    JoukowskyAirfoil(float(thickness), float(camber))
end

function JoukowskyAirfoil(thickness::Real, camber::Real)
    JoukowskyAirfoil(promote(thickness, camber)...)
end

leading_edge(airfoil::JoukowskyAirfoil) = position(airfoil, airfoil._LE_param)
trailing_edge(::JoukowskyAirfoil{T}) where {T} = Complex{T}(2)

function position(airfoil::JoukowskyAirfoil, t::Real)
    x = airfoil.μ + (1 - airfoil.μ) * cispi(2 * t)
    return x + inv(x)
end

is_clockwise(::JoukowskyAirfoil) = false


"""
    InterpolatedAirfoil{T} <: Airfoil{T}

An airfoil parameterization formed by interpolating a set of points.
"""
struct InterpolatedAirfoil{T} <: Airfoil{T}
    itp::AbstractInterpolation{Complex{T}}
    leading_edge::Complex{T}
    clockwise::Bool
    @doc """
        InterpolatedAirfoil(points, leading_edge, [clockwise])
    
    Interpolate the vector of complex `points` to form an airfoil with leading edge
    `leading_edge`. `clockwise` specifies whether `points` are ordered clockwsie, and is
    calculated if omitted.
    """
    function InterpolatedAirfoil(points::AbstractVector{Complex{T}},
                                 leading_edge::Complex{T};
                                 clockwise::Union{Bool, Nothing} = nothing
                                 ) where T <: AbstractFloat
        itp = LinearInterpolation(LinRange(0, 1, length(points)), points)
        if isnothing(clockwise)
            a = sum(zip(points, drop(points, 1))) do (z1, z2)
                (imag(z1) + imag(z2)) * (real(z2) - real(z1))
            end
            clockwise = a > 0
        end
        return new{T}(itp, leading_edge, clockwise)
    end
end

leading_edge(airfoil::InterpolatedAirfoil) = airfoil.leading_edge
trailing_edge(airfoil::InterpolatedAirfoil) = airfoil.itp(0)
position(airfoil::InterpolatedAirfoil, t::Real) = airfoil.itp(t)
is_clockwise(airfoil::InterpolatedAirfoil) = airfoil.clockwise


end # module
