def calculate_mass_loss_per_timestep(inputs):

    # set initial conditions
    n_layers = len(inputs["layer_thicknesses"])
    ks = k_in_subsurface_layers(inputs["k_star"], n_layers, inputs["layer_thicknesses"], inputs["extinction_coefficient"], inputs["delta_t"])
    qma = energy_at_surface(inputs["l_star"], inputs["qh"], inputs["qe"])
    masses = calculate_mass_per_layer(inputs["densities"], inputs["layer_thicknesses"])
    # set up output array
    total_mass_loss = []
    mass_loss_ablation = []
    mass_loss_internal = []

    for n in range(0,len(ks), 1):

        qmi = ks[n]
        p_ice = inputs["densities"][n]
        ma, mi, m = mass_loss_in_subsurface_layers(qma, qmi, p_ice, inputs["lf"])
        
        if n==0:
            print("here", m, mi, ma)
            total_mass_loss.append(m)
            mass_loss_ablation.append(ma)
            mass_loss_internal.append(mi)        
        else:
            total_mass_loss.append(mi)
            mass_loss_ablation.append(0)
            mass_loss_internal.append(mi)

    return total_mass_loss, mass_loss_ablation, mass_loss_internal


def calculate_mass_per_layer(densities, layer_thicknesses):
    masses = []
    for i in range(0, len(layer_thicknesses), 1):
        masses.append(densities[i] * layer_thicknesses[i])
    return masses

def k_in_subsurface_layers(k_star, n_layers, layer_thicknesses, extinction_coefficient, delta_t):
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

    for n in range(0, n_layers, 1):
        zl = sum(layer_thicknesses[0:n+1])
        zu = zl-layer_thicknesses[n]
        k_layer = ((k_star**(-extinction_coefficient*zu)) - (k_star**(-extinction_coefficient*zl))) * delta_t
        ks.append(k_layer)

    #print("Ks: \n", ks)
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
    # add ablation mass loss to upper layer only
    return ma, mi, m


def calculate_densities_in_each_layer(masses, total_mass_loss, mass_loss_internal):
    
    densities = (masses - total_mass_loss) / ()


    return densities


# def calculate_volume_lost_to_ablation(mass_loss_ablation, p_ice):

#     vas = []
#     for i in range(0, len(mass_loss_ablation), 1):
#         va = mass_loss_ablation[i] / p_ice[i]
#         vas.append(va)
#     return vas
