from util import *
from model import *
import csv
from config import inputs
import matplotlib.pyplot as plt

def run_model(inputs, config):
    
    # create instance of inputs that will update in each timestep
    inputs_i = inputs[0].copy()
    start = config["start"]
    stop = config["stop"]
    interval = config["interval"]
    csv_path = config["path_to_met_data"]
    outputs = []

    for i in range(start,stop,interval):
        
        if config["read_csv"]:
            inputs_i = load_inputs_from_csv(csv_path, inputs_i, i)

        inputs_i = set_initial_conditions(inputs_i)
        inputs_i = update_densities(inputs_i)
        inputs_i = replenish_lost_mass(inputs_i)
        inputs_i = calculate_density_at_t_plus_one(inputs_i)
        inputs_i = check_and_correct_negative_densities(inputs_i)

        # snapshot current state of inputs_i and push it to outputs
        outputs.append(inputs_i.copy())

        if config["outputs_to_csv"]:
            write_outputs_to_csv(outputs, config)


    return outputs

def write_outputs_to_csv(outputs, config):
    header = outputs[0].keys()
    
    with open(config["output_csv_savepath"], 'w') as f:
        data = []
        writer = csv.writer(f)
        writer.writerow(header)
        for i in range(0, len(outputs), 1):
            for key in outputs[i].keys():
                data.append(outputs[i][key])
            #data = [outputs[i]['day'],outputs[i]['time'],outputs[i]['qe'],outputs[i]['qh'],outputs[i]['densities'],outputs[i]['masses'], outputs[i]["volumes"]]
            writer.writerow(data)
            data = []
    return



def load_inputs_from_csv(filepath, inputs, timestep):

    with open(filepath) as csv_file:
        csv_file = csv.reader(csv_file, delimiter=',')
        rows = list(csv_file)
        inputs["day"] = float(rows[timestep+1][0])
        inputs["time"] = float(rows[timestep+1][1])
        inputs["inswrad"] = float(rows[timestep+1][2]) #timestep + 1 because row 0 is header
        inputs["avp"] = float(rows[timestep+1][3])
        inputs["airtemp"] = float(rows[timestep+1][4])
        inputs["windspd"] = float(rows[timestep+1][5])
        inputs["l_star"] = float(rows[timestep+1][6])
        inputs["rel_hum"] = float(rows[timestep+1][7])

    return inputs


def plot_density_profile(outputs, config):
    savepath = config["figure_savepath"]
    depths = outputs[0]["layer_thicknesses"]
    plt.figure(1)
    plt.ylim = (0.6,0.8)

    for i in range(0, len(outputs), 1):
        time = f"t{i}"
        plt.plot(depths, outputs[i]["densities"], label=time)

    plt.xlabel("depth below surface, m")
    plt.ylabel("density, kg m-3")
    plt.legend()
    plt.savefig(savepath)    
    plt.show()

    return
