from util import *
from model import *
import csv
from inputs import inputs
import matplotlib.pyplot as plt

def run_model(inputs, start, stop, interval, csv, csv_path):
    inputs_i = inputs[0].copy()
    outputs = []

    for i in range(start,stop,interval):
        inputs_i = inner(inputs_i, csv, csv_path, i)
        outputs.append(inputs_i)

    return outputs

def inner(inputs, csv, csv_path, i):
    
    inputs = inputs.copy()
    if csv:
        inputs = load_inputs_from_csv(csv_path, inputs, i)
    inputs = set_initial_conditions(inputs, i)
    inputs = update_densities(inputs, i)
    inputs = replenish_lost_mass(inputs, i)
    inputs = calculate_density_at_t_plus_one(inputs, i)
    inputs = setup_next_timestep(inputs, i)

    return inputs

def load_inputs_from_csv(filepath, inputs, timestep):

    with open(filepath) as csv_file:
        csv_file = csv.reader(csv_file, delimiter=',')
        rows = list(csv_file)
        inputs["k_star"] = float(rows[timestep+1][3]) #timestep + 1 because row 0 is header
        inputs["l_star"] = float(rows[timestep+1][4])
        inputs["wind_spd"] = float(rows[timestep+1][5])
        inputs["air_temp"] = float(rows[timestep+1][6])
        inputs["rel_hum"] = float(rows[timestep+1][7])

    return inputs


def plot_density_profile(outputs):
    depths = outputs[0]["layer_thicknesses"]
    plt.figure(1)
    plt.ylim = (0.6,0.8)

    for i in range(0, len(outputs), 1):
        time = f"t{i}"
        plt.plot(depths, outputs[i]["densities"], label=time)
        plt.legend()

    
    # for i in range(0, len(outputs), 1):
    #     plt.plot(depths, outputs[i]["densities"])
    
    plt.show()

    return
