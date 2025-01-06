


def run_model():
    for i in range():
        calculate_melt_per_timestep()
    return


def calculate_melt_per_timestep():
    # set initial conditions
    layer_thicknesses = [0.01, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.1] # meters
    n_layers = len(layer_thicknesses)
    extinction_coefficient = 0.006
    k_star = 100
    delta_t = 300
    densities = [0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750]
    lf = 333.7

    ks = k_per_layer(k_star, n_layers, layer_thicknesses, extinction_coefficient, delta_t)

    # set up output array
    ms = []
    for n in range(0,len(ks), 1):
        print(n)
        qma = k_star
        qmi = ks[n]
        p_ice = densities[n]
        ma, mi, m = mass_loss_per_layer(qma, qmi, p_ice, lf)
        ms.append(m)

    print("melt volumes: ", ms)
    
    return m


def k_per_layer(k_star, n_layers, layer_thicknesses, extinction_coefficient, delta_t):
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
        zu = sum(layer_thicknesses[0:n])-layer_thicknesses[n]
        zl = sum(layer_thicknesses[0:n])
        print("upper and lower layer boundaries", zu, zl)
        k_layer = (k_star**(extinction_coefficient*zu) - k_star**(extinction_coefficient*zl)) * delta_t
        ks.append(k_layer)
    return ks


def mass_loss_per_layer(qma, qmi, p_ice, lf):
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


run_model(1,1,1)
