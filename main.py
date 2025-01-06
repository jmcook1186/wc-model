
from util import *
from inputs import inputs

def run_model():
    for i in range(0,10,1):
        total_mass_loss, mass_loss_internal = calculate_mass_loss_per_timestep(inputs)
        masses = calculate_mass_per_layer(inputs)
        densities = calculate_densities_in_each_layer(masses, total_mass_loss, mass_loss_internal, inputs)
        inputs["densities"]= densities
        print(inputs["densities"])
    return

run_model()
