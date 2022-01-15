module PanelMethod

export Airfoils, AirfoilTools, Sheets, Sheet, SheetGroup, Flows, VortexStrengths, Flow

include("Utils.jl")
include("Airfoils.jl")
include("AirfoilTools.jl")

include("Sheets.jl")
include("Flows.jl")

import .Sheets
using .Sheets: Sheet, SheetGroup

import .Flows
using .Flows: VortexStrengths, Flow

end # module
