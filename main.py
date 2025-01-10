
from model import *
from util import *
from inputs import inputs


outputs = run_model(inputs, 0,10,1, True, '/home/joe/Desktop/MetEG(Sheet1).csv')

for i in range(0, len(outputs), 1):
    print(outputs[i]["densities"])

plot_density_profile(outputs)
