from util import *
from model import *
import csv
from inputs import inputs

def run_model(inputs, start, stop, interval, csv, csv_path):

    for i in range(start,stop,interval):
        if csv:
            inputs = load_inputs_from_csv(csv_path, inputs, i)
        
        print(inputs)
        set_initial_conditions(inputs, i)
        update_densities(inputs, i)
        replenish_lost_mass(inputs, i)
        calculate_density_at_t_plus_one(inputs, i)
        setup_next_timestep(inputs, i)

        # print(inputs[i])
        # print()

    return

def load_inputs_from_csv(filepath, inputs, timestep):

    with open(filepath) as csv_file:
        csv_file = csv.reader(csv_file, delimiter=',')
        rows = list(csv_file)
        inputs[timestep]["k_star"] = float(rows[timestep+1][3]) #timestep + 1 because row 0 is header
        inputs[timestep]["l_star"] = float(rows[timestep+1][4])
        inputs[timestep]["wind_spd"] = float(rows[timestep+1][5])
        inputs[timestep]["air_temp"] = float(rows[timestep+1][6])
        inputs[timestep]["rel_hum"] = float(rows[timestep+1][7])

    return inputs
