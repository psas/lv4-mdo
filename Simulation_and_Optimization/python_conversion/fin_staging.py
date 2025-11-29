
"""
Fin geometry is defined by the root, tip, sweep angle, semispan and thickness. We create a second rocket with new fin parameters to simulate staging. The weight reduction is also accounted for.
"""

# function to change fin parameters

def fin_resizing(LV4, root, tip, sweep, span, thickness, mass_red):
    new_rocket = LV4
    new_rocket.fin.root = root
    new_rocket.fin.tip = tip
    new_rocket.fin.sweep_angle = sweep
    new_rocket.fin.semispan = span
    new_rocket.fin.thickness = thickness
    new_rocket.mass = new_rocket.mass - mass_red
    return new_rocket



# given 0.025 time step calculate the time at which we should drop the fins
def find_drop_time(stability_list):
    closeness = []
    for item in stability_list:
        closeness.append(abs(item - 2.0))
    index = closeness.index(min(closeness))
    return (index * 0.025)


