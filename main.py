
from util import *
from inputs import inputs

def run_model():
    for i in range(0,10,1):
        m = calculate_mass_loss_per_timestep(inputs)
        print(m)
    return

run_model()
