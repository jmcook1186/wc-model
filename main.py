
from util import *
from inputs import inputs

def run_model():

    set_initial_conditions(inputs, 0)
    m, ma, mi = calculate_mass_of_melt(inputs, 0)
    print("MASSES: ", inputs[0]["masses"])
    print("mass loss: ", m)
    
    d = update_densities(inputs, 0)
    replenish_lost_mass(inputs, 0)
    print(inputs[0]["new_densities"])
    print(inputs[0]["new_masses"])

    # set_n_layers(inputs)
    # set_cumulative_depths(inputs)
    # set_ks(inputs)
    # set_masses(inputs)
    # set_volumes(inputs)

    # for i in range(0,3,1):
    #     total_mass_loss, mass_loss_ablation , mass_loss_internal = calculate_mass_loss_per_timestep(inputs)
    #     calculate_densities_in_each_layer(total_mass_loss, mass_loss_ablation, mass_loss_internal, inputs)
    
    return

run_model()
