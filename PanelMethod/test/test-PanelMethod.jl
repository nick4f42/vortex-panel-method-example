using Test
using PanelMethod

@testset verbose=true "PanelMethod" begin
    Test.detect_ambiguities(PanelMethod, recursive=true)

    include("test-Airfoils.jl")
end
