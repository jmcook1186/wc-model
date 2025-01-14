# wc-model

Python implementation a modified version of Schuster's (2001) weatherign crust development model.

## Quick start

Only one package is required that's not bundled in the Python standard library - matplotlib.
Install using `pip install matplotlib`. 

Then run the model from the project root directory using 

```
python main.py
```

This calls the `run_model()` function with the default config. You can change the config to alter how the model runs.

The default config has the `csv` parameter set to False, which means the model gets all its input data from the `inputs` dictionary defined in `inputs.py`. You can adjust the values in there as you see fit. This mode of operation should be sufficient if you are only doing single model runs or have coarse resolution data that you want to update manually.


## Met data

Often, you'll want to input data from a csv file. In this case, set `csv` to `True` int he call to `run_model()` and provide a filepath. Your csv file should have a header row and be structured as follows:

|TS	|DOY	|TOD	|k_star	|l_star	|wind_spd	|air_temp|	rel_hum|
|---|---|---|---|---|---|---|---|
|0	|217|	0	|10.85	|-56.475|	2.71|	0.73|	0.850775|
|1	|217|	100	|1.082	|-42.03	| 3.18	|0.2	|0.867466833|
|2	|217|	200	|1.322	|-39.269|	3.16|	0.17|	0.879960833|

...

## Modification from the Schuster model

We attempted to provide a faithful implementation of the Schuster (2001) model but the thesis where it is described lacks descriptions of some steps that I believe to be required for the model to yield density profiles in correct units.

Specifically, the equations 4.4 and 4.5, calculating Mi and Ma (mass loss due to internal weathering and mass loss due to ablation) are incomplete. In order to yield mass loss in kg, we replace the original equations with new ones.

The originals are

```
Mi = Qmi  / P * Lf

Ma = Qma / P * Lf
```

where Qmi = energy available for internal melting in Wm-2, Qma = energy available for surface melting in Wm-2, P = density of ice in kg m-3, Lf = latent heat of fusion for ice at 0 C (333700 J kg-1).

The updated equations are:

```
Ma = m * (Qma / m * Lf)
Mi = m * (Qmi / m * Lf)
```

Where m = mass of ice layer in kg, Qmi = energy available for internal melting in Wm-2, Qma = energy available for surface melting in Wm-2, , Lf = latent heat of fusion for ice at 0 C (333700 J kg-1)

Ma and Mi then represent the mass of ice lost to ablation.
