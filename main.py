
from model import *
from util import *
from inputs import inputs

# '/home/joe/Desktop/MetEG(Sheet1).csv'

outputs = run_model(inputs, 20,40,1, True, '/home/joe/Desktop/met_ilu22_hr.csv')

# for i in range(0, len(outputs), 1):
#     print(outputs[i]["densities"], "\n")

plot_density_profile(outputs)
