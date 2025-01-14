inputs = [{
    "layer_thicknesses" : [0.01, 0.02, 0.04, 0.08, 0.16, 0.32], # meters
    "area_w": 1, # m
    "area_l": 1, # m
    "dense_ice_constant": 890, # what density to use for underlying unweathered ice
    "extinction_coefficient" : 6, # m
    "k_star" : 380, # Wm-2
    "delta_t" : 3600, # seconds
    "densities" : [600, 600, 600, 600, 600, 600], # kgm-3
    "lf" : 333700, # J kg-1
    "l_star" : 100, # net long wave radiation
    "qh" : 60, # turbulent sensible heat flux density
    "qe" : 60 # turbulent latent heat flux density
}]
