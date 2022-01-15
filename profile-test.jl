using ProfileView

using PanelMethod
using .Airfoils

function dostuff()
    airfoil = JoukowskyAirfoil(0.1, 0.1) |> Xform(0, 1)

    leading_edge(airfoil)
    trailing_edge(airfoil)
    chordvec(airfoil)
    chordline(airfoil, 0.5)
    Airfoils.position(airfoil, 0.2)

    flow = Flow(Sheet(airfoil, 200), airfoil, 1, deg2rad(5))

    Sheets.aerodynamic_center(flow)
    Flows.lift_coeff(flow)
    Flows.center_of_pressure(flow)
    
    z = complex(3, 3)
    Flows.velocity(flow, z)
    Flows.streamfunction(flow, z)
end

dostuff()  # compile

ProfileView.@profview for _ in 1:50 dostuff() end
