using Test
using PanelMethod

@test isempty(Test.detect_ambiguities(PanelMethod, recursive=true))

include("test-Airfoils.jl")
