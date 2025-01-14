from util import *
from model import *
import csv
from inputs import inputs
import matplotlib.pyplot as plt

def run_model(inputs, start, stop, interval, csv, csv_path):
    
    # create instance of inputs that will update in each timestep
    inputs_i = inputs[0].copy()
    outputs = []

    for i in range(start,stop,interval):
        
        if csv:
            inputs_i = load_inputs_from_csv(csv_path, inputs_i, i)

        inputs_i = set_initial_conditions(inputs_i)
        inputs_i = update_densities(inputs_i)
        inputs_i = replenish_lost_mass(inputs_i)
        inputs_i = calculate_density_at_t_plus_one(inputs_i)
        inputs_i = check_and_correct_negative_densities(inputs_i)

        # snapshot current state of inputs_i and push it to outputs
        outputs.append(inputs_i.copy())

    return outputs


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

    
    # for i in range(0, len(outputs), 1):
    #     plt.plot(depths, outputs[i]["densities"])
    plt.xlabel("depth below surface, m")
    plt.ylabel("density, kg m-3")
    plt.legend()
    plt.savefig('/home/Desktop/wc-out.png')    
    plt.show()

    return
