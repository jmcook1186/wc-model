
from util import *
from inputs import inputs

def run_model():

    set_n_layers(inputs)
    set_ks(inputs)
    set_masses(inputs)
    set_volumes(inputs)

    for i in range(0,3,1):
        total_mass_loss, mass_loss_internal = calculate_mass_loss_per_timestep(inputs)
        densities = calculate_densities_in_each_layer(total_mass_loss, mass_loss_internal, inputs)
        print(densities)
    
    return

run_model()
