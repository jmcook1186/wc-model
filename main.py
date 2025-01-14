
from model import *
from util import *
from inputs import inputs

outputs = run_model(inputs, 1,10,1, False, '')

for i in range(0, len(outputs), 1):
    print(outputs[i]["densities"], "\n")

plot_density_profile(outputs)
