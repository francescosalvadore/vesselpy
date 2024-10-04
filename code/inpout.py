from stl_lib import read_stl, save_stl
from utils import decode_range_input, myround, print_log
import json
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0

def read_input():
    # instantiate
    config = ConfigParser()
    # parse existing file
    config.read('config.ini')
    # read values from general section
    gravity         = float(config.get('general', 'gravity'))
    p_atmo          = float(config.get('general', 'p_atmo'))
    refine_factor   = decode_range_input(config.get('general', 'refine_factor'))
    if config.has_option('general', 'max_trim'):
        max_trim        = float(config.get('general', 'max_trim'))
    else:
        max_trim = 20.
    if config.has_option('general', 'max_heel'):
        max_heel        = float(config.get('general', 'max_heel'))
    else:
        max_heel = 80.
    case            = config.get('general', 'case')
    # read values from case_* section
    case_section    = "case_"+case
    mode            = config.get(case_section, 'mode')
    filepath        = config.get(case_section, 'filepath')
    rho             = decode_range_input(config.get(case_section, 'rho'))
    mass            = decode_range_input(config.get(case_section, 'mass'))
    cog_x           = decode_range_input(config.get(case_section, 'cog_x'))
    cog_y           = decode_range_input(config.get(case_section, 'cog_y'))
    cog_z           = decode_range_input(config.get(case_section, 'cog_z'))
    z_water_fixed   = decode_range_input(config.get(case_section, 'z_water_fixed'))
    alfa_fixed      = decode_range_input(config.get(case_section, 'alfa_fixed'))
    heel_fixed      = decode_range_input(config.get(case_section, 'heel_fixed'))
    cor_x           = config.get(case_section, 'cor_x')
    cor_x           = ["cog"] if cor_x == "cog" else [float(cor_x)] 
    cor_y           = config.get(case_section, 'cor_y')
    cor_y           = ["cog"] if cor_y == "cog" else [float(cor_y)] 
    cor_z           = config.get(case_section, 'cor_z')
    cor_z           = ["cog"] if cor_z == "cog" else [float(cor_z)] 
    heel_stability  = decode_range_input(config.get(case_section, 'heel_stability'))

    #cog_z           = [cog_z[0]]  # cog_z role to be investigated 

    stl_obj = read_stl("stl/"+filepath)
    # The choice of 0 for blocked quantites should be improved. Especially cog_x for instance.
    # Note: for unknown DoFs an iteration is needed and then "auto" allows to start the initial
    # iteration value in a reasonable manner. Otherways, explicit values act as starting values.
    # For known DoFs the values act as real values and auto cannot be used clearly.
    if   mode == "weight_lcg_tcg":    # 3DOF
        if z_water_fixed == "auto": 
            z_water_fixed = [0.5*(stl_obj["boundingbox"][2] + stl_obj["boundingbox"][5])]
        if alfa_fixed == "auto":
            alfa_fixed    = [0]
        if heel_fixed == "auto":
            heel_fixed    = [0]
    elif mode == "weight_lcg_heel":   # 2DOF
        if z_water_fixed == "auto": 
            z_water_fixed = [0.5*(stl_obj["boundingbox"][2] + stl_obj["boundingbox"][5])]
        if alfa_fixed == "auto": 
            alfa_fixed    = [0]
        if cog_y == "auto": # useless, there is no iteration on cog, it is given or found in one shot
            cog_y         = [0.5*(stl_obj["boundingbox"][1] + stl_obj["boundingbox"][4])]
    elif mode == "weight_trim_tcg":   # 2DOF
        if z_water_fixed == "auto": 
            z_water_fixed = [0.5*(stl_obj["boundingbox"][2] + stl_obj["boundingbox"][5])]
        if cog_x == "auto": # useless
            cog_x         = [0.5*(stl_obj["boundingbox"][0] + stl_obj["boundingbox"][3])]
        if heel_fixed == "auto": 
            heel_fixed    = [0] 
    elif mode == "weight_trim_heel":  # 1DOF
        if z_water_fixed == "auto": 
            z_water_fixed = [0.5*(stl_obj["boundingbox"][2] + stl_obj["boundingbox"][5])]
        if cog_x == "auto": # useless
            cog_x         = [0.5*(stl_obj["boundingbox"][0] + stl_obj["boundingbox"][3])]
        if cog_y == "auto": # useless
            cog_y         = [0.5*(stl_obj["boundingbox"][1] + stl_obj["boundingbox"][4])]
    elif mode == "sinkage_trim_heel": # 0DOF
        # Beware: even if z_water, trim and heel are given, CoR plays a crucial role to
        # evaluate position and z_water interpretation as a consequence.
        mass              = ["auto"]
        if cog_x == "auto": # useless
            cog_x         = [0.5*(stl_obj["boundingbox"][0] + stl_obj["boundingbox"][3])]
        if cog_y == "auto": # useless
            cog_y         = [0.5*(stl_obj["boundingbox"][1] + stl_obj["boundingbox"][4])]

    if config.has_option(case_section, 'control_points'):
        control_points = json.loads(config.get(case_section,'control_points'))
    else:
        control_points = []

    return mode, gravity, p_atmo, refine_factor, max_trim, max_heel, case, filepath, rho, \
           mass, cog_x, cog_y, cog_z, z_water_fixed, alfa_fixed, heel_fixed, \
           cor_x, cor_y, cor_z, control_points, heel_stability

def save_results(results, folder):
    fr = open(folder+"/results.dat","w")
    for k,v in results.items():
        fr.write(f"{k} = {myround(v,3)}\n")
    fr.close()

def save_water_stl(stl_obj, z_water_opt, folder):
    xmin, ymin, zmin, xmax, ymax, zmax = stl_obj["boundingbox"]
    xmin_t = xmin-(xmax-xmin)*0.8 
    ymin_t = ymin-(ymax-ymin)*0.8
    xmax_t = xmax+(xmax-xmin)*0.8 
    ymax_t = ymax+(ymax-ymin)*0.8
    pminmin = [xmin_t, ymin_t, z_water_opt]
    pminmax = [xmin_t, ymax_t, z_water_opt]
    pmaxmin = [xmax_t, ymin_t, z_water_opt]
    pmaxmax = [xmax_t, ymax_t, z_water_opt]
    t1 = [pminmax, pminmin, pmaxmin]
    t2 = [pminmax, pmaxmax, pmaxmin]
    stl_obj_water = {}
    stl_obj_water["triangles"] = [t1, t2]
    stl_obj_water["normals"] = [[0,0,1], [0,0,1]]
    stl_obj_water["centers"] = [[0,0,0], [0,0,0]] # not used
    save_stl(stl_obj_water, folder+"/z_water_opt.stl")
