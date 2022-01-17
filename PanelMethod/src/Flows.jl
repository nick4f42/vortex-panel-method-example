module Flows

using Base.Iterators: drop

using ..Utils
using ..Airfoils
using ..Sheets
using ..Sheets: _lengths, _point_count

export
    VortexStrengths, Flow, freestream, angle_of_attack, velocity, streamfunction,
    circulation, pressure_coeff, ndim_force, lift_coeff, moment_coeff, center_of_pressure


const SubVector = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true} where T


"""
    VortexStrengths(sheetgroup::SheetGroup)
    
Storage for the vortex strengths necessary to describe `sheetgroup`.
"""
struct VortexStrengths{T<:AbstractFloat}
    flat::Vector{T}
    vec::Vector{SubVector{T}}
    function VortexStrengths(sheetgroup::SheetGroup{T}) where T
        γ_flat = Vector{T}(undef, _point_count(sheetgroup))

        next_i = 1
        γ_vec = Vector{SubVector{T}}(undef, length(sheetgroup))
        for (k, sheet) in zip(eachindex(γ_vec), sheetgroup)
            i = next_i
            next_i += _point_count(sheet)
            γ_vec[k] = view(γ_flat, i : next_i-1)
        end

        return new{T}(γ_flat, γ_vec)
    end
end


"""
    Flow(sheets, chord::AbstractChord, [γ::VortexStrengths,] v::Number)
    
Create a flow from `sheets` with complex freestream velocity `v` and store the resulting
vortex strengths in `γ`. `chord` determines how quantities are nondimensionalized. `sheets`
may be a single [`Sheet`](@ref), vector of [`Sheet`](@ref), or [`SheetGroup`](@ref).

    Flow(sheets, chord::AbstractChord, [γ::VortexStrengths,] v::Real, α::Real)

Create a flow with freestream velocity of magnitude `v` and angle of attaack `α`.
"""
struct Flow{T, G<:SheetGroup{T}, C<:AbstractChord}
    sheetgroup::G
    chord::C
    γ::VortexStrengths{T}
    u_inf::Complex{T}
    function Flow(sheetgroup::G, chord::C, γ::VortexStrengths{T}, v::Number
                  ) where {T, G<:SheetGroup{T}, C<:AbstractChord}
        u_inf::Complex{T} = v
        @. γ.flat = real(u_inf) * sheetgroup._x_re + imag(u_inf) * sheetgroup._x_im
        return new{T, G, C}(sheetgroup, chord, γ, u_inf)
    end
end

function Flow(sheetgroup::SheetGroup{T}, chord::AbstractChord, γ::VortexStrengths{T},
              v::Real, α::Real) where T
    return Flow(sheetgroup, chord, γ, v * cis(α) * chordunit(chord))
end

function Flow(sheetgroup::SheetGroup, chord::AbstractChord, args::Vararg{Number, N}
              ) where N
    return Flow(sheetgroup, chord, VortexStrengths(sheetgroup), args...)
end

function Flow(sheets, chord::AbstractChord, args::Vararg{Any, N}) where N
    Flow(SheetGroup(sheets), chord, args...)
end

"""
    freestream(flow::Flow)

The complex freestream velocity of `flow`.
"""
freestream(flow::Flow) = flow.u_inf

"""
    angle_of_attack(flow::Flow)

The angle of attack of `flow`.
"""
angle_of_attack(flow::Flow) = angle(conj(chordvec(flow.chord)) * freestream(flow))

"""
    velocity(flow::Flow, z)

The flow velocity at `z`.
"""
function velocity(flow::Flow{T}, z::Complex{T}) where T
    u::Complex{T} = 0

    for (s, γ) in zip(flow.sheetgroup, flow.γ.vec)
        Δu::Complex{T} = 0
        for (z1, z2, γ1, γ2) in zip(s.endpoints, drop(s.endpoints, 1), γ, drop(γ, 1))
            tmp1 = z1 - z2
            tmp2 = 1im * abs(tmp1) / (2 * pi * tmp1^2)
            tmp3 = log((z - z1) / (z - z2)) * tmp2
            tmp4 = tmp1 * tmp2
        
            c1 = tmp3 * (z2 - z) - tmp4
            c2 = tmp3 * (z - z1) + tmp4
            Δu += conj(c1 * γ1 + c2 * γ2)
        end
        u += Δu
    end

    return u + freestream(flow)
end

function velocity(flow::Flow{T}, z::Number) where T
    velocity(flow, convert(Complex{T}, z))
end

"""
    streamfunction(flow::Flow, z)

The stream function at `z`.
"""
function streamfunction(flow::Flow{T}, z::Complex{T}) where T
    ψ::T = 0
    for (s, γ) in zip(flow.sheetgroup, flow.γ.vec)
        Δψ::T = 0
        for (z1, z2, γ1, γ2) in zip(s.endpoints, drop(s.endpoints, 1), γ, drop(γ, 1))
            tmp1 = z - z1
            tmp2 = z - z2
            tmp3 = z1 - z2
            tmp4 = γ1 * (2 * z2 - z1 - z) + γ2 * tmp1

            Δψ += abs(tmp3) / 8pi * (
                2 * (real(tmp4 / tmp3) + real(tmp1 * tmp4 / tmp3^2 * log(tmp1 / tmp2)))
                + γ1 - γ2 + (γ1 + γ2) * log(abs2(tmp2))
            )
        end
        ψ += Δψ
    end

    return ψ + ccross(freestream(flow), z)
end

function streamfunction(flow::Flow{T}, z::Number) where T
    streamfunction(flow, Complex{T}(z))
end

"""
    circulation(flow::Flow)

The clockwise circulation around all the sheets in `flow`.
"""
function circulation(flow::Flow{T}) where T
    Γ::T = 0
    for (sheet, γ) in zip(flow.sheetgroup, flow.γ.vec)
        ΔΓ::T = 0
        for (L, γ1, γ2) in zip(_lengths(sheet), γ, drop(γ, 1))
            ΔΓ += L * (γ1 + γ2) / 2
        end
        Γ += ΔΓ
    end
    return Γ
end

"""
    pressure_coeff(flow::Flow, z)

The pressure coefficient at `z`.
"""
pressure_coeff(flow::Flow, z::Number) = 1 - abs2(velocity(flow, z)) / abs2(freestream(flow))

"""
    panel_ndim_force(flow::Flow)

The total complex pressure force on every sheet found by summing that of each panel.
Nondimensionalized by chordlength.
"""
function panel_ndim_force(flow::Flow{T}) where T
    F::Complex{T} = 0
    for sheet in flow.sheetgroup
        ΔF::Complex{T} = 0
        for (p, n, L) in zip(sheet.ctrl_points, sheet.ctrl_normals, _lengths(sheet))
            ΔF += -n * pressure_coeff(flow, p) * L
        end
        F += ΔF
    end
    return F / chordlength(flow.chord)  # nondimensionalize
end

"""
    ndim_force(flow::Flow)

The total complex force found by assuming the Kutta-Joukowski theorem.
"""
ndim_force(flow::Flow) = 1im * freestream(flow) * lift_coeff(flow)

"""
    lift_coeff(flow::Flow)

The total lift coefficient found by assuming the Kutta-Joukowski theorem.
"""
function lift_coeff(flow::Flow)
    return 2 * circulation(flow) / (abs(freestream(flow)) * chordlength(flow.chord))
end

"""
    lift_coeff(flow)

The moment coefficient about `z` found by summing pressure forces.
"""
function moment_coeff(flow::Flow{T}, z::Complex{T}) where T
    M::T = 0
    for sheet in flow.sheetgroup
        ΔM::T = 0
        for (p, n, L) in zip(sheet.ctrl_points, sheet.ctrl_normals, _lengths(sheet))
            ΔM += ccross(p - z, -n * pressure_coeff(flow, p) * L)
        end
        M += ΔM
    end
    # clockwise-positive moment convention
    return -M / chordlength(flow.chord)^2  # nondimensionalize
end

function moment_coeff(flow::Flow{T}, z::Number) where T
    moment_coeff(flow, convert(Complex{T}, z))
end

"""
    center_of_pressure(flow::Flow)

Return the center of pressure as a fraction of chord length of `flow.chord`.
"""
function center_of_pressure(flow::Flow)
    C_m = moment_coeff(flow, leading_edge(flow.chord))
    return -C_m / ccross(chordvec(flow.chord), ndim_force(flow))
end

"""
    aerodynamic_center(flow::Flow)

The aerodynamic center of `flow`'s underlying `SheetGroup` and chord.
"""
function Sheets.aerodynamic_center(flow::Flow)
    Sheets.aerodynamic_center(flow.sheetgroup, flow.chord)
end

end # module
