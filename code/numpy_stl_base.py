import numpy as np
from stl import mesh

def get_features_from_of(filename, mass, cog, logfile_name="progress.log"):
    # Using an existing closed stl file:
    your_mesh = mesh.Mesh.from_file(filename)
    
    volume, cog_computed, inertia_tensor = your_mesh.get_mass_properties()
    density = mass / volume
    #print("Volume                                           = {0}".format(volume))
    #print("Position of the center of gravity (COG)          = {0}".format(cog))
    #print("Position of the center of gravity (COG_computed) = {0}".format(cog_computed))
    #print("Inertia matrix at expressed at the COG_computed  = {0}".format(inertia_tensor[0,:]))
    #print("                                                   {0}".format(inertia_tensor[1,:]))
    #print("                                                   {0}".format(inertia_tensor[2,:]))
 
    inertia_tensor = density*inertia_tensor
    inertia_tensor = [float(x) for x in inertia_tensor.ravel()]

    return float(volume), [float(x) for x in cog_computed], inertia_tensor

if __name__ == "__main__":
    volume, cog_computed, inertia_tensor = get_features_from_of("../../../RUN/DTMB/stl/DTMB_V1.stl",500.,[2.5,-0.005,0.3])
    print("input: ",500, [2.5,-0.005,0.3])
    print(volume, cog_computed, inertia_tensor)
