module AirfoilsTest

using Test

using PanelMethod
using .Airfoils
using .Airfoils: position

macro test_chord(chord)
    quote
        chord = $(esc(chord))

        @inferred leading_edge(chord)
        @inferred trailing_edge(chord)
        @inferred chordvec(chord)

        @inferred chordunit(chord)
        @test abs(chordunit(chord)) ≈ 1

        @inferred chordlength(chord)

        @inferred chordline(chord, 0.5)
        @inferred chordline(chord, 0)
        @inferred chordline(chord, 1//1)
    end
end

macro test_airfoil(airfoil)
    quote
        airfoil = $(esc(airfoil))

        @test_chord airfoil

        chordlen = chordlength(airfoil)
        
        @inferred position(airfoil, 0)
        @inferred position(airfoil, 0.5)
        @inferred position(airfoil, 1//4)

        @test position(airfoil, 0) ≈ position(airfoil, 1) atol=1e-8*chordlen
    end
end

@testset "Airfoils" begin

    @testset "Xform and mirror" begin
        @test_chord @inferred Xform(0, 1)
        @test_chord @inferred Xform(1//2, 3.0 + 4.0im)
        @test_chord @inferred Xform(1 - 3im, 4 + 5im)
        @test_chord @inferred Xform(JoukowskyAirfoil(0.1, 0.1))

        @test_airfoil @inferred JoukowskyAirfoil(0.1, 0.1) |> Xform(0, 1)
        @test_airfoil @inferred JoukowskyAirfoil(0.1, 0.1) |> Xform(1//2, 1//2)

        @test_airfoil @inferred JoukowskyAirfoil(0.1, 0.1) |> mirror
        @test_airfoil @inferred JoukowskyAirfoil(0.1, 0.1) |> Xform(0, 1) |> mirror
    end

    @testset "airfoil types" begin
        @testset "JoukowskyAirfoil" begin
            @test_airfoil @inferred JoukowskyAirfoil(0.1, 0.1)
            @test_airfoil @inferred JoukowskyAirfoil(1//8, 1//4)
            @test_airfoil @inferred JoukowskyAirfoil(0.1, 1//4)

            @test_throws DomainError JoukowskyAirfoil(-0.1, 0.2)
            @test_throws DomainError JoukowskyAirfoil(1//8, -0.1)
        end

        @testset "InterpolatedAirfoil" begin
            let points = position.(Ref(JoukowskyAirfoil(0.1, 0.1)), LinRange(0, 1, 10))
                @test_airfoil @inferred InterpolatedAirfoil(points, points[cld(end, 2)])
            end
            let points = ComplexF64[1, 1im, -1, -1im, 1]
                @test_airfoil @inferred InterpolatedAirfoil(points, points[cld(end, 2)])
            end
        end
    end

end

end # module
