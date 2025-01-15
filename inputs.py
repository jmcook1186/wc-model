inputs = [{
    # variables
    "day": 250,
    "time": 100,
    "avp": 660,
    "airtemp": 5.612,
    "windspd": 3.531,
    "lat": 0,
    "lon": 0,
    "lon_ref" : 0,
    "summertime" : 0,
    "slope" : 1.,
    "aspect" : 90.,
    "elevation" : 1000.,
    "albedo" : 0.35,
    "roughness" : 0.005,
    "met_elevation" : 1000.,
    "lapse" : 0.0065,
    "layer_thicknesses" : [0.01, 0.02, 0.04, 0.08, 0.16, 0.32], # meters
    "area_w": 1, # m
    "area_l": 1, # m
    "dense_ice_constant": 890, # what density to use for underlying unweathered ice
    "extinction_coefficient" : 3, # m
    "inswrad" : 380, # Wm-2
    "delta_t" : 3600, # seconds
    "densities" : [600, 600, 600, 600, 600, 600], # kgm-3
    "lf" : 333700, # J kg-1
    "l_star" : 100, # net long wave radiation

}]
