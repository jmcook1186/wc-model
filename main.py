
from util import *
from inputs import inputs

def run_model():

    set_n_layers(inputs)
    set_cumulative_depths(inputs)
    set_ks(inputs)
    set_masses(inputs)
    set_volumes(inputs)


    for i in range(0,1,1):
        total_mass_loss, mass_loss_ablation , mass_loss_internal = calculate_mass_loss_per_timestep(inputs)
        calculate_densities_in_each_layer(total_mass_loss, mass_loss_ablation, mass_loss_internal, inputs)
    

    return

run_model()
