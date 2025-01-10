from math import exp
import matplotlib.pyplot as plt


def setup_next_timestep(inputs, timestep):
    inputs["densities"] = inputs["next_densities"]
    inputs["masses"] = inputs["new_masses"]
    return inputs

def set_initial_conditions(inputs, timestep):
    inputs["n_layers"] = set_n_layers(inputs, timestep)
    inputs["volumes"] = set_volumes(inputs, timestep)
    inputs["masses"] = set_mass(inputs, timestep)
    inputs["total_depth"] = set_total_column_depth(inputs, timestep)
    inputs["k"] = set_k(inputs, timestep)
    return inputs

def set_mass(inputs, timestep):
    masses = []
    for i in range(0, len(inputs["densities"]), 1):
        masses.append(inputs["densities"][i] * inputs["volumes"][i])
    return masses

def set_n_layers(inputs, timestep):
    return len(inputs["layer_thicknesses"])

def set_volumes(inputs, timestep):
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


def set_total_column_depth(inputs,timestep):
    return sum(inputs["layer_thicknesses"])


def set_k(inputs, timestep):
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
        k_star_n = (k_star*(exp(-extinction_coefficient * zu))-(k_star*(exp(-extinction_coefficient*zl))))*delta_t
        k_stars.append(k_star_n)
    
    # account for condition when K* is positive but qma is negative (pg 114, between eqs 4.12 and 4.13)
    qma = calculate_surface_melt_energy(inputs, timestep)
    if (k_stars[0] > 0) & (qma < 0):
        k_stars[0] = k_star + qma
    
    return k_stars

def calculate_surface_melt_energy(inputs, timestep):
    return inputs["l_star"] + inputs["qh"] + inputs["qe"]

def calculate_mass_of_melt(inputs, timestep):
    qma = calculate_surface_melt_energy(inputs, timestep)
    ma = qma / (inputs["densities"][0] * inputs["lf"])
    mi=[]
    
    for i in range(0, len(inputs["k"]), 1):
        mi.append(inputs["k"][i] / (inputs["densities"][i] * inputs["lf"]))
    
    m = mi
    m[0] = ma + mi[0]

    return m, ma, mi


def update_densities(inputs, timestep):

    densities_before = inputs["densities"]
    masses = inputs["masses"]
    volumes = inputs["volumes"]
    densities_after = []

    
    m, ma, mi = calculate_mass_of_melt(inputs, timestep)
    
    # top layer only
    densities_after.append((masses[0] - m[0])/(volumes[0]-(ma/densities_before[0])))
    
    # subsurface layers
    for i in range(1, len(densities_before),1):
        densities_after.append((masses[i] - mi[i])/(volumes[i]))

    inputs["new_densities"] = densities_after
    return inputs


def replenish_lost_mass(inputs, timestep):
    new_densities = inputs["new_densities"]
    volumes = inputs["volumes"]
    m, ma, mi = calculate_mass_of_melt(inputs, timestep)
    new_masses = []
    
    for i in range(0, len(new_densities), 1):
        if ma > 0:  # if there is ablation
            if i < len(new_densities)-1:
                new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (new_densities[i+1]*(ma/new_densities[i])))
            else:
                new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (0.89*(ma/new_densities[i])))
        else: # if there's no ablation, just give current value
            new_masses.append(new_densities[i]*(volumes[i]))

    inputs["new_masses"] = new_masses
    
    return inputs

def calculate_density_at_t_plus_one(inputs, timestep):
    new_densities = inputs["new_densities"]
    volumes = inputs["volumes"]
    m, ma, mi = calculate_mass_of_melt(inputs, timestep)

    next_densities =[]
    
    for i in range(0, len(new_densities)-1, 1):
        if new_densities[i] > 0: #as long as there's no complete collapse
            next_densities.append((new_densities[i]*(1-(ma/new_densities[i]/volumes[i])))+ (new_densities[i+1]*(ma/new_densities[i]/volumes[i])))
        else:
            next_densities.append(new_densities[i+1])

    # handle bottom layer
    if new_densities[-1] > 0:
        next_densities.append((new_densities[i]*(1-(ma/new_densities[i]/volumes[i])))+ (0.89*(ma/new_densities[i]/volumes[i])))
    else:
        next_densities.append(0.89)

    inputs["next_densities"] = next_densities

    return inputs
