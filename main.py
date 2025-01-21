
from model import *
from util import *
from config import inputs, config

outputs = run_model(inputs, config)

# for i in range(0, len(outputs), 1):
#     print(i, ": ", outputs[i]["qh"], "\n")

plot_density_profile(outputs, config)
