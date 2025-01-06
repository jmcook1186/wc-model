
from util import *
from inputs import inputs

def run_model():
    for i in range(0,1,1):
        m, ma, mi = calculate_mass_loss_per_timestep(inputs)
    
    print(m, "\n", mi, "\n", ma)
    return

run_model()
