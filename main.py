
from model import *
from util import *
from inputs import inputs

outputs = run_model(inputs, 1,30,1, True, '/home/joe/Desktop/Samplmet.csv')

for i in range(0, len(outputs), 1):
    print(i, ": ", outputs[i]["qh"], "\n")

plot_density_profile(outputs)
