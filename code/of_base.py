#!/usr/bin/env python3

from utils import myrun, verify_openfoam

def get_features_from_of(filename, mass, cog, logfile_name="progress.log"):
    # verify openfoam installation
    verify_openfoam()
    
    out, err = myrun('surfaceInertia '+filename+' | grep "Mass:" ')
    volume   = float(out.split()[1])
    density  = mass/volume
    
    cog      = '('+str(cog[0])+' '+str(cog[1])+' '+str(cog[2])+')'
    #out, err = myrun('surfaceInertia -referencePoint \''+cog+'\' -density '+str(density)+' '+filename+' | grep "Centre of mass:" ')
    out, err = myrun('surfaceInertia -referencePoint \''+cog+'\' -density '+str(density)+' '+filename+' >> output.of_base')
    out, err = myrun('cat output.of_base | grep "Centre of mass:"')
    a = out.split(":")[1]
    #print(out)
    b = a.strip("()\n ")
    cog_computed = b.split()
    cog_computed = [float(x) for x in cog_computed]

    out, err = myrun('cat output.of_base | grep -A 1 "Inertia tensor around centre of mass:" | tail -n 2')
    #print(out)
    a = out.split(":")[1]
    b = a.strip("()\n ")
    inertia_tensor = b.split()
    inertia_tensor = [float(x) for x in inertia_tensor]

    out, err = myrun('cat output.of_base >> '+logfile_name)
    out, err = myrun('rm output.of_base')
    out, err = myrun('mv axes.obj work')

    # RIMETTERE SENZA
    #cog_computed_str = '('+str(cog_computed[0])+' '+str(cog_computed[1])+' '+str(cog_computed[2])+')'
    #print("aaa: ",cog_computed_str)
    #out, err = myrun('surfaceInertia -referencePoint \''+cog_computed_str+'\' -density '+str(density)+' '+filename+' >> output.of_base')
    #out, err = myrun('cat output.of_base')
    #print(out)
    # RIMETTERE SENZA
    return volume, cog_computed, inertia_tensor

if __name__ == "__main__":
    #volume, cog_computed, inertia_tensor = get_features_from_of("stl/w3_real.stl",32.6,[0.5695,0.,0.156])
    volume, cog_computed, inertia_tensor = get_features_from_of("../../../RUN/DTMB/stl/DTMB_V1.stl",500.,[2.5,-0.005,0.3])
    print(volume, cog_computed, inertia_tensor)
