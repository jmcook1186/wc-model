

def calculate_mass_loss_per_timestep(inputs):

    # set initial conditions
    n_layers = len(inputs["layer_thicknesses"])
    ks = k_in_subsurface_layers(inputs["k_star"], n_layers, inputs["layer_thicknesses"], inputs["extinction_coefficient"], inputs["delta_t"])
    qma = energy_at_surface(inputs["l_star"], inputs["qh"], inputs["qe"])
    # set up output array
    ms = []

    for n in range(0,len(ks), 1):
        qmi = ks[n]
        p_ice = inputs["densities"][n]
        ma, mi, m = mass_loss_in_subsurface_layers(qma, qmi, p_ice, inputs["lf"])
        ms.append(m)

    print("melt mass per layer (kg): ", ms)
    
    return m


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
        k_layer = (k_star**(extinction_coefficient*zu) - k_star**(extinction_coefficient*zl)) * delta_t
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
    m = ma + mi
    return ma, mi, m
