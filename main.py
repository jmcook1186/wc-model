
from util import *

def run_model():
    for i in range(0,10,1):
        m = calculate_melt_per_timestep()
        print(m)
    return

run_model()
