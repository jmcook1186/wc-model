from math import exp
import matplotlib.pyplot as plt

def set_n_layers(inputs):
    inputs["n_layers"] = len(inputs["layer_thicknesses"])
    print("setting n-layers:")
    print(inputs["n_layers"])
    print()
    return

def set_cumulative_depths(inputs):
    cumulative_depths = []
    for n in range(0, len(inputs["layer_thicknesses"]), 1):
        print(n)
        if n == 0:
            cumulative_depths.append(inputs["layer_thicknesses"][0])
        else:
            cumulative_depths.append(sum(inputs["layer_thicknesses"][0:n]) + inputs["layer_thicknesses"][n])
    
    inputs["cumulative_depths"] = cumulative_depths

    print("setting cumulative depths")
    print(cumulative_depths)
    print()
    return

def set_ks(inputs):
    print("K* = ", inputs["k_star"])
    ks = k_in_subsurface_layers(inputs)
    inputs["ks"] = ks
    print("K*n = ", ks)
    return 

def calculate_mass_loss_per_timestep(inputs):

    qma = energy_at_surface(inputs["l_star"], inputs["qh"], inputs["qe"])
    print("Qma = ", qma)
    # set up output array
    total_mass_loss = []
    mass_loss_ablation = []
    mass_loss_internal = []

    for n in range(0,len(inputs["ks"]), 1):

        qmi = inputs["ks"][n]
        p_ice = inputs["densities"][n]
        ma, mi, m = mass_loss_in_subsurface_layers(qma, qmi, p_ice, inputs["lf"])
        print(m, mi, ma)

        if n==0:
            total_mass_loss.append(m)
            mass_loss_ablation.append(ma)
            mass_loss_internal.append(mi)        
        else:
            total_mass_loss.append(mi)
            mass_loss_ablation.append(0)
            mass_loss_internal.append(mi)

    return total_mass_loss, mass_loss_ablation, mass_loss_internal

def set_masses(inputs):
    masses = []
    for i in range(0, len(inputs["layer_thicknesses"]), 1):
        masses.append(inputs["densities"][i] * inputs["layer_thicknesses"][i])
    inputs["masses"] = masses
    return

def k_in_subsurface_layers(inputs):
    """K* value for each layer in column using Beer's law

    Parameters:
    k_star: surface radiative flux in Wm-2
    n_layers: number of layers in column
    layer_thicknesses: array of layer thicknesses in m
    extinction_coefficient: extinction coefficient at each layer in m

    Returns:
    ks: array of fluxes available for internal melting (Qmi)

    """
    ks = []
    k_star = inputs["k_star"]
    n_layers = inputs["n_layers"]
    layer_thicknesses = inputs["layer_thicknesses"]
    extinction_coefficient = inputs["extinction_coefficient"]
    delta_t = inputs["delta_t"]

    for n in range(0, n_layers, 1):
        zl = sum(layer_thicknesses[0:n+1])
        zu = zl-layer_thicknesses[n]
        k_layer = ((k_star * exp(-extinction_coefficient*zu)) - (k_star*exp(-extinction_coefficient*zl))) * delta_t
        ks.append(k_layer)
    return ks

def energy_at_surface(l_star, qh, qe):
    """calculates mass loss at upper surface layer

    Parameters:
    l_star: net longwave radiation Wm-2
    qh: turbulent sensible heat flux density Wm-2
    qe: turbulent latent heat flux density Wm-2

    Returns:
    qma: mass loss in surface layer
    """
    return l_star + qh + qe

def mass_loss_in_subsurface_layers(qma, qmi, p_ice, lf):
    """Mass of ice lost per horizontal layer

    Parameters:
    qma: surface lowering flux (Wm-2)
    qmi: internal melting flux (Wm-2)
    p_ice = density of the ice layer (kg m-3)
    lf = latenty heat of fusion (333.7 j g-1 at 0deg C)

    Returns:
    ma: mass lost due to surface lowering (kg)
    mi: mass lost due to internal melting (kg)
    m: total mass loss (kg) as sum of ma and mi

    """
    
    mi = qmi / (p_ice * lf)
    ma = qma /(p_ice * lf)
    m = mi + ma

    return ma, mi, m

def calculate_densities_in_each_layer(total_mass_loss, mass_loss_ablation, mass_loss_internal, inputs):
    densities = []
    masses = []
    volumes = []

    for n in range(0,len(inputs["densities"]),1):
        
        p = inputs["densities"][n]
        volume = inputs["masses"][n] / p
        volumes.append(volume)
        volume_loss_ablation = total_mass_loss[n] / p

        if n == 0:
            temp_density = (inputs["masses"][n] - total_mass_loss[n]) / (volume - volume_loss_ablation)
            masses.append(temp_density*(inputs["volumes"][n] - volume_loss_ablation)+(inputs["densities"][n]*volume_loss_ablation))
            density = p*(1-(volume_loss_ablation/volume)+(temp_density*(volume_loss_ablation/volume)))
            if density < 0:
                densities.append(inputs["densities"][n+1])
            else:
                densities.append(density)
            
        if n < len(inputs["densities"])-1:
            print("n", n)
            temp_density =(inputs["masses"][n] - mass_loss_internal[n]) / volume
            masses.append(temp_density*(inputs["volumes"][n] - volume_loss_ablation)+(p*volume_loss_ablation))
            density = p*(1-(volume_loss_ablation/volume)+(temp_density*(volume_loss_ablation/volume)))
            
            if density < 0:
                densities.append(inputs["densities"][n+1])
            else:
                densities.append(density)

        else:
            if density < 0:
                densities.append(0.89)
            else:
                densities.append(density)

    # update initial conditions
    inputs["densities"] = densities
    inputs["masses"] = masses
    inputs["volumes"] = volumes

    print("DENSITIES: ", densities)
    return 

def set_volumes(inputs):
    volumes = []
    for i in range(0, len(inputs["layer_thicknesses"]), 1):
        volumes.append(inputs["masses"][i] / inputs["densities"][i])
    inputs["volumes"] = volumes
    return


def replenish_ice_from_below(inputs):

    return
