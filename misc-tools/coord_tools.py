"Tools to convert Corsika coordinates to GROPTICS coordinates and center them around the array center. Useful for get_telescope_X and get_telescope_Y functions."

import numpy as np

#=============== 3-TEL ARRAY================#

COORDINATES_CORSIKA_3 = {
    "Heli": (-22.20e2, -76.58e2, 39.70e2),
    "Fern": (97.56e2, 11.55e2, 34.66e2),
    "Winter": (-75.36e2, 67.04e2, 49.17e2),
} # in cm

COORDINATES_GROPTICS_3 = {
    name: (-y, x, z) for name, (x, y, z) in COORDINATES_CORSIKA_3.items()
}

coords_array_3 = np.array(list(COORDINATES_GROPTICS_3.values()), dtype=float)
center_3 = coords_array_3.mean(axis=0)

COORDINATES_CENTERED_3 = {
    name: tuple(coord - center_3) for name, coord in COORDINATES_GROPTICS_3.items()
}

def get_telescope_X_3(n): # in cm
    if n==1:
        return COORDINATES_CENTERED_3["Heli"][0]/100
    elif n==2:
        return COORDINATES_CENTERED_3["Fern"][0]/100
    else: #n==3
        return COORDINATES_CENTERED_3["Winter"][0]/100  

def get_telescope_Y_3(n): # in cm
    if n==1:
        return COORDINATES_CENTERED_3["Heli"][1]/100
    elif n==2:
        return COORDINATES_CENTERED_3["Fern"][1]/100
    else: #n==3
        return COORDINATES_CENTERED_3["Winter"][1]/100

#=============== 7-TEL ARRAY================#

COORDINATES_CORSIKA_7 = {
    "Heli":    (-0.10e2,   97.18e2,   39.70e2),
    "Winter":  (-220.15e2, 33.67e2,   49.17e2),
    "Fern":    (-130.43e2, 190.40e2,  34.66e2),
    "Gattini": (177.58e2,  -333.33e2, 42.82e2),
    "Tower":   (-402.48e2, 133.47e2,  45.23e2),#
    "Antler":  (132.97e2,  -39.01e2,  37.22e2),
    "Vent":    (0.10e2,    -97.18e2,  46.35e2),
} # in cm

COORDINATES_GROPTICS_7 = {
    name: (-y, x, z) for name, (x, y, z) in COORDINATES_CORSIKA_7.items()
}

coords_array_7 = np.array(list(COORDINATES_GROPTICS_7.values()), dtype=float)
center_7 = coords_array_7.mean(axis=0)

COORDINATES_CENTERED_7 = {
    name: tuple(coord - center_7) for name, coord in COORDINATES_GROPTICS_7.items()
}

def get_telescope_X_7(n): # in meters
    if n==1:
        return COORDINATES_CENTERED_7["Heli"][0]/100
    elif n==2:
        return COORDINATES_CENTERED_7["Winter"][0]/100
    elif n==3:
        return COORDINATES_CENTERED_7["Fern"][0]/100
    elif n==4:
        return COORDINATES_CENTERED_7["Gattini"][0]/100
    elif n==5:
        return COORDINATES_CENTERED_7["Tower"][0]/100
    elif n==6:
        return COORDINATES_CENTERED_7["Antler"][0]/100
    else: #n==7
        return COORDINATES_CENTERED_7["Vent"][0]/100  

def get_telescope_Y_7(n): # in meters
    if n==1:
        return COORDINATES_CENTERED_7["Heli"][1]/100
    elif n==2:
        return COORDINATES_CENTERED_7["Winter"][1]/100
    elif n==3:
        return COORDINATES_CENTERED_7["Fern"][1]/100
    elif n==4:
        return COORDINATES_CENTERED_7["Gattini"][1]/100
    elif n==5:
        return COORDINATES_CENTERED_7["Tower"][1]/100
    elif n==6:
        return COORDINATES_CENTERED_7["Antler"][1]/100
    else: #n==7
        return COORDINATES_CENTERED_7["Vent"][1]/100