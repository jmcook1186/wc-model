from math import exp
import matplotlib.pyplot as plt


def set_initial_conditions(inputs, timestep):
    inputs[timestep]["n_layers"] = set_n_layers(inputs, timestep)
    inputs[timestep]["volumes"] = set_volumes(inputs, timestep)
    inputs[timestep]["masses"] = set_mass(inputs, timestep)
    inputs[timestep]["total_depth"] = set_total_column_depth(inputs, timestep)
    inputs[timestep]["k"] = set_k(inputs, timestep)
    return

def set_mass(inputs, timestep):
    masses = []
    for i in range(0, len(inputs[timestep]["densities"]), 1):
        masses.append(inputs[timestep]["densities"][i] * inputs[timestep]["volumes"][i])
    return masses

def set_n_layers(inputs, timestep):
    return len(inputs[timestep]["layer_thicknesses"])

def set_volumes(inputs, timestep):
    volumes =[]
    for t in inputs[timestep]["layer_thicknesses"]:
        volumes.append(t*inputs[timestep]["area_l"]*inputs[timestep]["area_w"])
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
    return sum(inputs[timestep]["layer_thicknesses"])


def set_k(inputs, timestep):
    """K* value for each layer in column using Beer's law

    Parameters:
    inputs: array of dicts with input values for each timestep

    Returns:
    ks: array of fluxes available for internal melting (Qmi)

    """

    k_star = inputs[timestep]["k_star"]
    layer_thicknesses = inputs[timestep]["layer_thicknesses"]
    extinction_coefficient = inputs[timestep]["extinction_coefficient"]
    delta_t = inputs[timestep]["delta_t"]

    k_stars = []

    for i in range(0, len(layer_thicknesses), 1):
        zl = sum(layer_thicknesses[0:i+1])
        zu = zl-layer_thicknesses[i]
        k_star_n = (k_star*(exp(-extinction_coefficient * zu))-(k_star*(exp(-extinction_coefficient*zl))))*delta_t
        k_stars.append(k_star_n)
    
    # account for condition when K* is positive but qma is negative (pg 114, between eqs 4.12 and 4.13)
    qma = calculate_surface_melt_energy(inputs, timestep)
    if k_stars[0] > 0 & qma < 0:
        k_stars[0] = k_star + qma
    
    return k_stars

def calculate_surface_melt_energy(inputs, timestep):
    return inputs[timestep]["l_star"] + inputs[timestep]["qh"] + inputs[timestep]["qe"]

def calculate_mass_of_melt(inputs, timestep):
    qma = calculate_surface_melt_energy(inputs, timestep)
    ma = qma / (inputs[timestep]["densities"][timestep] * inputs[timestep]["lf"])
    mi=[]
    for i in range(0, len(inputs[timestep]["k"]), 1):
        mi.append(inputs[timestep]["k"][i] / (inputs[timestep]["densities"][i] * inputs[timestep]["lf"]))
    
    m = mi
    m[0] = ma + mi[0]

    return m, ma, mi

def update_densities(inputs, timestep):
    densities_before = inputs[timestep]["densities"]
    masses = inputs[timestep]["masses"]
    volumes = inputs[timestep]["volumes"]
    densities_after = []
    m, ma, mi = calculate_mass_of_melt(inputs, timestep)
    
    # top layer only
    densities_after.append((masses[0] - m[0])/(volumes[0]-(ma/densities_before[0])))
    
    # subsurface layers
    for i in range(1, len(densities_before),1):
        densities_after.append((masses[i] - mi[i])/(volumes[i]))

    inputs[timestep]["new_densities"] = densities_after
    return densities_after


def replenish_lost_mass(inputs, timestep):
    new_densities = inputs[timestep]["new_densities"]
    volumes = inputs[timestep]["volumes"]
    m, ma, mi = calculate_mass_of_melt(inputs, timestep)
    new_masses = []
    for i in range(0, len(new_densities), 1):
        if i < len(new_densities)-1:
            new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (new_densities[i+1]*(ma/new_densities[i])))
        else:
            new_masses.append(new_densities[i]*(volumes[i]-(ma/new_densities[i])) + (0.89*(ma/new_densities[i])))
    inputs[timestep]["new_masses"] = new_masses
    
    return

def calculate_density_at_t_plus_one():
    
    return
# def calculate_mass_loss_per_timestep(inputs):

#     qma = energy_at_surface(inputs["l_star"], inputs["qh"], inputs["qe"])
#     # set up output array
#     total_mass_loss = []
#     mass_loss_ablation = []
#     mass_loss_internal = []

#     for n in range(0,len(inputs["ks"]), 1):

#         qmi = inputs["ks"][n]
#         p_ice = inputs["densities"][n]

#         ma, mi, m = mass_loss_in_subsurface_layers(qma, qmi, p_ice, inputs["lf"])

#         if n==0:
#             total_mass_loss.append(m)
#             mass_loss_ablation.append(ma)
#             mass_loss_internal.append(mi)        
#         else:
#             total_mass_loss.append(mi)
#             mass_loss_ablation.append(0)
#             mass_loss_internal.append(mi)

#     return total_mass_loss, mass_loss_ablation, mass_loss_internal

# def set_masses(inputs):
#     masses = []
#     for i in range(0, len(inputs["layer_thicknesses"]), 1):
#         masses.append(inputs["densities"][i] * inputs["layer_thicknesses"][i])
#     inputs["masses"] = masses
#     return

# def k_in_subsurface_layers(inputs):
    # """K* value for each layer in column using Beer's law

    # Parameters:
    # k_star: surface radiative flux in Wm-2
    # n_layers: number of layers in column
    # layer_thicknesses: array of layer thicknesses in m
    # extinction_coefficient: extinction coefficient at each layer in m

    # Returns:
    # ks: array of fluxes available for internal melting (Qmi)

    # """
    # ks = []
    # k_star = inputs["k_star"]
    # n_layers = inputs["n_layers"]
    # layer_thicknesses = inputs["layer_thicknesses"]
    # extinction_coefficient = inputs["extinction_coefficient"]
    # delta_t = inputs["delta_t"]

    # for n in range(0, n_layers, 1):
    #     zl = sum(layer_thicknesses[0:n+1])
    #     zu = zl-layer_thicknesses[n]
    #     k_layer = ((k_star * exp(-extinction_coefficient*zu)) - (k_star*exp(-extinction_coefficient*zl))) * delta_t
    #     ks.append(k_layer)

    # qma = energy_at_surface(inputs["l_star"], inputs["qh"], inputs["qe"])
    # print(qma, k_star+qma, k_star)
    
    # # account for condition when K* is positive but qma is negative (pg 114, between eqs 4.12 and 4.13)
    # if qma <0:
    #     ks[0] = (((k_star + qma) * exp(-extinction_coefficient*zu)) - ((k_star + qma)*exp(-extinction_coefficient*zl))) * delta_t
    
    # return ks

# def energy_at_surface(l_star, qh, qe):
#     """calculates mass loss at upper surface layer

#     Parameters:
#     l_star: net longwave radiation Wm-2
#     qh: turbulent sensible heat flux density Wm-2
#     qe: turbulent latent heat flux density Wm-2

#     Returns:
#     qma: mass loss in surface layer
#     """
#     return l_star + qh + qe

# def mass_loss_in_subsurface_layers(qma, qmi, p_ice, lf):
#     """Mass of ice lost per horizontal layer

#     Parameters:
#     qma: surface lowering flux (Wm-2)
#     qmi: internal melting flux (Wm-2)
#     p_ice = density of the ice layer (kg m-3)
#     lf = latenty heat of fusion (333.7 j g-1 at 0deg C)

#     Returns:
#     ma: mass lost due to surface lowering (kg)
#     mi: mass lost due to internal melting (kg)
#     m: total mass loss (kg) as sum of ma and mi

#     """
    
#     mi = qmi / (p_ice * lf)
#     ma = qma /(p_ice * lf)
#     m = mi + ma

#     return ma, mi, m


# def calculate_densities_in_each_layer(total_mass_loss, mass_loss_ablation, mass_loss_internal, inputs):

#     densities = []
#     masses = []
#     volumes = []
#     next_densities = []
#     next_volumes =  []
#     next_masses = []

#     # density calculation in top layer
#     top_layer_density = (inputs["masses"][0] - total_mass_loss[0])/((inputs["masses"][0]/inputs["densities"][0])-(total_mass_loss[0]/inputs["densities"][0]))
#     densities.append(top_layer_density)

#     for n in range(0,len(inputs["densities"])-1,1):

#         density_before = inputs["densities"][n]
#         volume_before = inputs["masses"][n]/density_before
#         mass_before = inputs["masses"][n]

#         # if we are in bottom layer
#         if n == len(inputs["densities"])-1:
#             density_after = (mass_before - mass_loss_internal / volume_before)
#             mass_after = density_after*(volume_before - (mass_loss_ablation/density_before))
#             volume_after = mass_after/density_after
#             if density_after < 0:
#                 densities.append(0.89)
#             else:
#                 densities.append(density_after)
#             masses.append(mass_after)
#             densities.append(density_after)
#             volumes.append(volume_after)

#         else:
#             density_after = (mass_before - mass_loss_internal[n] / volume_before)
#             mass_after = density_after*(volume_before - (mass_loss_ablation[n]/density_before))
#             volume_after = mass_after/density_after
            
#             if density_after < 0:
#                 densities.append(inputs["densities"][n+1])
#             else:
#                 densities.append(density_after)

#             masses.append(mass_after)
#             densities.append(density_after)
#             volumes.append(volume_after)


#     # for n in range(0,len(inputs["densities"])-2,1):
#     #     if qma > 0:
#     #         density_after = densities[n]
#     #         volume_after = volumes[n]
#     #         mass_after = masses[n]

#     #         new_mass = density_after*(volume_before-(mass_loss_ablation[n]/density_before))+(densities[n+1]*(mass_loss_ablation[n]/density_before))
#     #         new_density = new_mass/volume_after
#     #         next_densities.append(new_density)
#     #         next_masses.append(new_mass)

#     # # deal with bottom layer again
#     # next_densities.append(density_after)
#     # next_masses.append(mass_after)
#     # next_volumes.append(volume_after)

#     # update initial conditions
#     # inputs["densities"] = next_densities
#     # inputs["masses"] = next_masses
#     # inputs["volumes"] = next_volumes

#     inputs["densities"] = densities
#     inputs["masses"] = masses
#     inputs["volumes"] = volumes

#     print("DENSITIES: ", densities)

#     return 

# def set_volumes(inputs):
#     volumes = []
#     for i in range(0, len(inputs["layer_thicknesses"]), 1):
#         volumes.append(inputs["masses"][i] / inputs["densities"][i])
#     inputs["volumes"] = volumes
#     return


# def replenish_ice_from_below(inputs):

#     return
