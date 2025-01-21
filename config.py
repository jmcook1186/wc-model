inputs = [{
    # variables
    "day": 100,
    "time": 1000,
    "windspd": 1,
    "rel_hum": 10,
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
    "extinction_coefficient" : 6, # m
    "delta_t" : 3600, # seconds
    "densities" : [600, 600, 600, 600, 600, 600], # kgm-3
    "lf" : 333700, # J kg-1
    "l_star" : 100, # net long wave radiation
}]

config = {
    "read_csv": True,
    "output_fluxes_to_csv": True,
    "flux_csv_savepath": "/home/joe/Desktop/flux_data.csv",
    "path_to_met_data": '/home/joe/Desktop/Samplmet.csv',
    "start": 1,
    "stop": 20,
    "interval": 1,
    "figure_savepath": "/home/joe/Desktop/wc-model.png"
}
