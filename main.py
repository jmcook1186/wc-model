
from util import *
from inputs import inputs

def run_model():

    for i in range(0,2,1):
        set_initial_conditions(inputs, i)
        update_densities(inputs, i)
        replenish_lost_mass(inputs, i)
        calculate_density_at_t_plus_one(inputs, i)
        setup_next_timestep(inputs, i)

        print(inputs[i])
        print()

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
