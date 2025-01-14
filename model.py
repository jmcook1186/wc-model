from math import exp
import matplotlib.pyplot as plt

def set_initial_conditions(inputs):
    inputs["n_layers"] = set_n_layers(inputs)
    inputs["volumes"] = set_volumes(inputs)
    inputs["masses"] = set_mass(inputs)
    inputs["total_depth"] = set_total_column_depth(inputs)
    inputs["k"] = set_k(inputs)

    return inputs

def set_mass(inputs):
    masses = []
    for i in range(0, len(inputs["densities"]), 1):
        masses.append(inputs["densities"][i] * inputs["volumes"][i])
    return masses

def set_n_layers(inputs):
    return len(inputs["layer_thicknesses"])

def set_volumes(inputs):
    volumes =[]
    for t in inputs["layer_thicknesses"]:
        volumes.append(t*inputs["area_l"]*inputs["area_w"])
    return volumes

def set_cumulative_depths(inputs, timestep):
    cumulative_depths = []
    for n in range(0, len(inputs["layer_thicknesses"]), 1):
        if n == 0:
            cumulative_depths.append(inputs["layer_thicknesses"][timestep])
        else:
            cumulative_depths.append(sum(inputs["layer_thicknesses"][0:n]) + inputs["layer_thicknesses"][n])
    return cumulative_depths

def set_total_column_depth(inputs):
    return sum(inputs["layer_thicknesses"])

def set_k(inputs):
    """K* value for each layer in column using Beer's law

    Parameters:
    inputs: array of dicts with input values for each timestep

    Returns:
    ks: array of fluxes available for internal melting (Qmi)

    """

    k_star = inputs["k_star"]
    layer_thicknesses = inputs["layer_thicknesses"]
    extinction_coefficient = inputs["extinction_coefficient"]
    delta_t = inputs["delta_t"]
    k_stars = []

    for i in range(0, len(layer_thicknesses), 1):
        zl = sum(layer_thicknesses[0:i+1])
        zu = zl-layer_thicknesses[i]
        ku = k_star*exp(-extinction_coefficient*zu)
        kl = k_star*exp(-extinction_coefficient*zl)
        kdiff = ku-kl
        k_star_n = kdiff * delta_t
        k_stars.append(k_star_n)
    
    # account for condition when K* is positive but qma is negative (pg 114, between eqs 4.12 and 4.13)
    qma = calculate_surface_melt_energy(inputs)
    if (k_stars[0] > 0) & (qma < 0):
        k_stars[0] = k_star + qma
    
    return k_stars

def calculate_surface_melt_energy(inputs):
    return inputs["l_star"] + inputs["qh"] + inputs["qe"]


def calculate_mass_of_melt(inputs):

    qma = calculate_surface_melt_energy(inputs)

    ma = inputs["masses"][0] * (qma / (inputs["masses"][0] * inputs["lf"]))

    if ma < 0:
        ma = 0
    
    mi=[]

    for i in range(0, len(inputs["k"]), 1):
        mi.append(inputs["masses"][i] * (inputs["k"][i] / (inputs["masses"][i] * inputs["lf"])))

    m = mi
    
    m[0] = ma + mi[0]

    return m, ma, mi


def update_densities(inputs):

    densities_before = inputs["densities"]
    masses = inputs["masses"]
    volumes = inputs["volumes"]

    densities_after = []
    m, ma, mi = calculate_mass_of_melt(inputs)

    # top layer only
    densities_after.append((masses[0] - m[0])/(volumes[0]-(ma/densities_before[0]))) # eq 4.6a
    
    # subsurface layers
    for i in range(1, len(densities_before),1):
        densities_after.append((masses[i] - mi[i])/(volumes[i])) # eq 4.6b

    inputs["densities"] = densities_after

    return inputs


def replenish_lost_mass(inputs):
    new_densities = inputs["densities"]
    volumes = inputs["volumes"]
    m, ma, mi = calculate_mass_of_melt(inputs)
    new_masses = []

    for i in range(0, len(new_densities), 1):
        if i ==0:
            if new_densities[i] <=0:
                # upper layer replenishment
                new_masses.append(new_densities[i]*(volumes[i]-0.01) + (0.01*new_densities[i+1]))
            else:
                new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (new_densities[i+1]*(ma/new_densities[i])))
        elif i < len(new_densities)-1:
            #middle-layers
            new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (new_densities[i+1]*(ma/new_densities[i])))
        else:
            # bottom layer
            new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (inputs["dense_ice_constant"]*(ma/new_densities[i])))

    inputs["masses"] = new_masses
    
    return inputs

def calculate_density_at_t_plus_one(inputs):

    volumes = inputs["volumes"]
    masses = inputs["masses"]

    next_densities =[]

    for i in range(0, len(masses), 1):
        next_densities.append(masses[i]/volumes[i])

    inputs["densities"] = next_densities
    return inputs


def check_and_correct_negative_densities(inputs):
    
    densities = inputs["densities"]
    
    for i in range(len(densities)-1, -1, -1):
        if densities[i] < 0:
            if i == len(densities)-1:
                densities[i] = inputs["dense_ice_constant"]
            else:
                densities[i] = densities[i-1] 

    inputs["densities"] = densities

    return inputs
