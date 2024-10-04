import copy
import os
import shutil
import itertools 

from stl_lib import read_stl, refine_stl, rotate_stl_y, rotate_stl_x, save_stl
from utils import print_log, myround
from inpout import read_input, save_water_stl, save_results

#from of_base import get_features_from_of
from numpy_stl_base import get_features_from_of

from pv_base import extract_wl, extract_immersed_hull
from solver import hydrostatic_solver, full_rotate_point, rotate_point, get_features_base 
from solver import compute_floatation, compute_variations, get_form_coefficients, compute_forces_moment
from solver import get_intact_stability_heel, get_intact_stability_trim, get_metacenter_heel_general, get_metacenter_trim_general

if __name__ == "__main__":

    # Input imported from config.ini
    mode, gravity, p_atmo, refine_factor_l, max_trim, max_heel, case, filepath, rho_l, \
        mass_l, cog_x_l, cog_y_l, cog_z_l, z_water_fixed_l, alfa_fixed_l, heel_fixed_l, \
        cor_x_l, cor_y_l, cor_z_l, control_points, heel_stability_l = read_input()

    if not os.path.exists("output"): os.mkdir("output")
    if not os.path.exists("work"): os.mkdir("work")

    # Loop over all cases
    i_sim = 0

    all_cases = list(itertools.product(refine_factor_l, rho_l, mass_l, cog_x_l, cog_y_l, cog_z_l, z_water_fixed_l, alfa_fixed_l, heel_fixed_l, cor_x_l, cor_y_l, cor_z_l))
    n_hydro_tot = len(all_cases)
    print_log("Total number of hydrostatics analyses to be performed: ",n_hydro_tot)

    for refine_factor, rho, mass, cog_x, cog_y, cog_z, z_water_fixed, alfa_fixed, heel_fixed, cor_x, cor_y, cor_z in all_cases:

        i_sim += 1
        print_log("===> Running simulation :",i_sim)
        folder = "output/"+str(i_sim)  
        if not os.path.exists(folder): os.mkdir(folder)

        if cor_x == "cog": cor_x = cog_x
        if cor_y == "cog": cor_y = cog_y
        if cor_z == "cog": cor_z = cog_z
        cog       = [cog_x, cog_y, cog_z]
        cor       = [cor_x, cor_y, cor_z]

        with open(folder+"/info.dat","w") as fi:
            for v in ["filepath", "max_trim", "max_heel", "mode", "rho", "mass", "cog_x", "cog_y", "cog_z", \
                      "z_water_fixed", "alfa_fixed", "heel_fixed", \
                      "cor", "cor_x", "cor_y", "cor_z", "heel_stability_l", "control_points"]:
                fi.write(f"{v} = {myround(globals()[v])}\n")

        results = dict()
        results["cor"] = cor

        # Hydrostatics
        print_log("[1/5] Computing hydrostatic equilibrium...", end="")
        stl_obj = read_stl("stl/"+filepath)
        if refine_factor > 1:
            refine_stl(stl_obj, refine_factor)
            save_stl(stl_obj, "work/hull_refined.stl")
        else:
            #print_log("Refine skipped")
            pass
        stl_obj_save = copy.deepcopy(stl_obj)
        z_water_opt, alfa_opt, heel_opt, mass_r, cog_x_r, cog_y_r, cog_z_r, equilibrium_found = \
            hydrostatic_solver(mode, stl_obj, mass, rho, p_atmo, gravity, max_trim, max_heel, cog, cor, z_water_fixed, alfa_fixed, heel_fixed)
        print_log("...completed")

        results["z_water_opt"], results["alfa_opt"], results["heel_opt"] = z_water_opt, alfa_opt, heel_opt
        results["mass_r"], results["cog_x_r"], results["cog_y_r"], results["cog_z_r"] = mass_r, cog_x_r, cog_y_r, cog_z_r
        results["displacement"] = mass_r*gravity

        # Position and save hydrostatic model. Beware: CoR is used as center of rotations during the
        # convergence of the equations. CoR can be different from CoG and surely it is when dealing 
        # with a mode with fixed DOFs since CoG is not known.
        stl_obj_rotate_x  = rotate_stl_x(stl_obj_save, heel_opt, cor)
        stl_obj_static    = rotate_stl_y(stl_obj_rotate_x, alfa_opt, cor)
        save_stl(stl_obj_static, "work/hull_rotated.stl")
        shutil.copyfile("work/hull_rotated.stl",folder+"/hull_rotated.stl")

        # save waterline stl
        save_water_stl(stl_obj_static, z_water_opt, folder)

        # Redefine mass, and cog-x/y according to fixed DoFs. 
        mass   = mass_r
        cog[0] = cog_x_r
        cog[1] = cog_y_r
        cog[2] = cog_z_r

        # Find the CoG in the reference frame of the initial geometry. Actually, CoR (center of rotation) can be 
        # different from CoG and surely it is for modes where some DoFs are fixed since CoG is unknown.
        # The found CoG refers to the hydrostatic reference frame and it has to be counter-rotated around 
        # CoR to find cog_ini_ref, i.e., CoG in the initial reference frame (initial geometry).
        # Compute z_water_opt_rot_cog considering imposed rotations around cog_ini_ref instead of cog_input
        # This quantity is useful for external usage of data. If I want to specify angles and z_water, with 
        # rotations done around cog while here rotations were done around provisional cog_input.
        print_log("cog recomputed: ",str(cog))
        cog_ini_ref                    = full_rotate_point(cor, cog, -alfa_opt, -heel_opt, first_rotation="trim")
        print_log("cog_ini_ref: ",str(cog_ini_ref))
        results["cog_ini_ref"]         = cog_ini_ref
        results["z_water_opt_rot_cog"] = z_water_opt - (cog[2]-cog_ini_ref[2])

        # Get position of control points in hydrostatic equilibrium (i.e., following rotations around CoR)
        control_points_r = []
        for icp,cp in enumerate(control_points):
            r_ty, r_tz = rotate_point((cor[1],cor[2]), (cp[1],cp[2]), heel_opt)
            r_tx, r_tz = rotate_point((cor[0],cor[2]), (cp[0],r_tz),  alfa_opt)
            cp_el = [r_tx, r_ty, r_tz]
            control_points_r.append(cp_el)
            results["control_points_"+str(icp)] = cp_el
        results["control_points_r"] = control_points_r

        # Get inertia tensor around cog
        print_log("[2/5] Computing inertia tensor and basic features...", end="")
        volume, cog_computed, inertia_tensor = get_features_from_of("work/hull_rotated.stl", mass, cog)
        results["total_volume"], results["cog_computed"], results["inertia_tensor"] = volume, cog_computed, inertia_tensor
        results["inertia_tensor_diag"] = [inertia_tensor[0], inertia_tensor[4], inertia_tensor[8]]
    
        # Get basic features
        LOA, BOA, D, T, xmin, ymin, zmin, xmax, ymax, zmax = get_features_base(stl_obj_save, stl_obj_static, z_water_opt)
        results["LOA"], results["BOA"], results["D"], results["T"] = LOA, BOA, D, T
        results["bounding_box_x"] = [xmin, xmax]
        results["bounding_box_y"] = [ymin, ymax]
        results["bounding_box_z"] = [zmin, zmax]
        print_log("...completed")
    
        # Floatation
        print_log("[3/5] Computing floatation, buoyancy, and form coefficients...", end="")
        extract_wl("work/hull_rotated.stl", "work/hull_wl.stl", z_water_opt)
        stl_wl_obj = read_stl("work/hull_wl.stl")
        AWP, LCF, TCF, LWL, BWL, wl_xmin, wl_xmax, wl_ymin, wl_ymax, mom_in_x, mom_in_y = compute_floatation(stl_wl_obj, z_water_opt)
        cof = [LCF, TCF, z_water_opt]
        results["cof"] = cof
        results["cof_mom_x"] = mom_in_x
        results["cof_mom_y"] = mom_in_y
        # LCF and TCF to be verified as reference coordinates
        cof_ini_ref = full_rotate_point(cor, cof, -alfa_opt, -heel_opt, first_rotation="trim")
        results["cof_ini_ref"] = cof_ini_ref
        #USELESSwl_min = [wl_xmin, wl_ymin, z_water_opt]
        #USELESSwl_max = [wl_xmax, wl_ymax, z_water_opt]
        #USELESSresults["wl_min_ini_ref"] = full_rotate_point(cor, wl_min, -alfa_opt, -heel_opt, first_rotation="trim")
        #USELESSresults["wl_max_ini_ref"] = full_rotate_point(cor, wl_max, -alfa_opt, -heel_opt, first_rotation="trim")
        #USELESSresults["wl_xmin"], results["wl_xmax"], results["wl_ymin"], results["wl_ymax"] = wl_xmin, wl_xmax, wl_ymin, wl_ymax

        results["AWP"], results["LCF"], results["TCF"], results["LWL"], results["BWL"] = AWP, cof_ini_ref[0], cof_ini_ref[1], LWL, BWL
    
        # Variations
        weight_to_sink, moment_to_trim = compute_variations(stl_obj_static, z_water_opt, cog, mass, LWL, LOA, rho, p_atmo, gravity)
        results["weight_to_sink"], results["moment_to_trim"] = weight_to_sink, moment_to_trim
    
        # Buoyancy
        extract_immersed_hull("work/hull_rotated.stl", "work/hull_immersed.stl", z_water_opt)
        displaced_volume, cob, _ = get_features_from_of("work/hull_immersed.stl", mass, cog)
        print_log("extract hull immersed, cob : ",z_water_opt, cob)
        results["displaced_volume"], results["cob"] = displaced_volume, cob
        cob_ini_ref = full_rotate_point(cor, cob, -alfa_opt, -heel_opt, first_rotation="trim")
        results["cob_ini_ref"] = cob_ini_ref
        results["LCB"], results["TCB"], results["VCB"] = cob_ini_ref[0], cob_ini_ref[1], cob_ini_ref[2]

        extract_immersed_hull("work/hull_rotated.stl", "work/hull_immersed_nowl.stl", z_water_opt, include_wl=False)
        stl_obj_immersed_nowl = read_stl("work/hull_immersed_nowl.stl")
        wet_area = sum(stl_obj_immersed_nowl["areas"])
        results["wet_area"] = wet_area
    
        # Form coefficients
        Cb, Cvp, Cwp, Cws = get_form_coefficients(LWL, BWL, T, displaced_volume, AWP, wet_area)
        results["Cb"], results["Cvp"], results["Cwp"], results["Cws"] = Cb, Cvp, Cwp, Cws
        print_log("...completed")
    
        # Static stability
        # basic method: BMt = I_y/V_immersed ; BMl = I_x/V_immersed 
        print_log("[4/5] Computing intact stability...", end="")
        results["BMt"] = mom_in_y/displaced_volume
        results["BMl"] = mom_in_x/displaced_volume
        results["GMt"] = results["BMt"] - cog[2] + cob[2]
        results["GMl"] = results["BMl"] - cog[2] + cob[2]
        results["Mt"] = results["GMt"]  + cog[2] - z_water_opt
        results["Ml"] = results["GMl"]  + cog[2] - z_water_opt
        results["metacenter_t"] = [cog[0], cog[1], cog[2]+results["GMt"]]
        results["metacenter_l"] = [cog[0], cog[1], cog[2]+results["GMl"]]

        # alternative method
        metacenter_rot = False
        if metacenter_rot:
            GMt, BMt, Mt, metacenter_t = get_intact_stability_heel(stl_obj_static, z_water_opt, cog, cob, cof, rho, p_atmo, gravity)
            results["GMt_rot"], results["BMt_rot"], results["Mt_rot"], results["metacenter_t_rot"] = GMt, BMt, Mt, metacenter_t
            GMl, BMl, Ml, metacenter_l = get_intact_stability_trim(stl_obj_static, z_water_opt, cog, cob, cof, rho, p_atmo, gravity)
            results["GMl_rot"], results["BMl_rot"], results["Ml_rot"], results["metacenter_l_rot"] = GMl, BMl, Ml, metacenter_l
    
        # general method
        metacenter_check = False
        if metacenter_check:
            metacenter = get_metacenter_heel_general(stl_obj_static, z_water_opt, mass, cog)
            results["metacenter_t_general_method"] = metacenter
            metacenter = get_metacenter_trim_general(stl_obj_static, z_water_opt, mass, cog)
            results["metacenter_l_general_method"] = metacenter
        print_log("...completed")

        # Save results
        save_results(results, folder)
    
        # Stability curve
        # Note that the 0 degree condition is taken as the original model orientation, not the
        # equilibrium flotation plane. If a non-zero TCG or a non-zero Model Heel are entered, there
        # will be a non-zero righting arm at 0 degrees of heel. Zero righting arm will correspond to the
        # heel angle at the equilibrium flotation plane.
        # Since the heel stability works using the initial stl geometry, cog_ini_ref must be used here.
        # This is different from cog_r if CoR is not CoG, which is surely true at least for modes with limited 
        # DoFs, because in those modes, CoG is not known and rotations are performed around CoR which is not CoG.
        # This means that the  hydrostatic position and initial condition do not share the cog.
        print_log("[5/5] Computing stability curve...",end="")
        fs = open(folder+"/stab.dat", "w")
        fs.write(f"# heel z_water_opt alfa_opt fz my mx by\n")
        for heel in heel_stability_l:
            print_log(f"heel={heel}...", end="")
            print(f"heel={heel}...", end="")
            # Since at this stage we for sure know where cog_ini_ref is, we can use it also as center of rotations to avoid drifting of cog.
            # It seems to be more reasonable to rotate around it instead of the given input cor.
            cor_stab = cog_ini_ref # cor

            z_water_opt, alfa_opt, heel_opt, mass_r, cog_x_r, cog_y_r, cog_z_r, equilibrium_found = \
                hydrostatic_solver("weight_lcg_heel", stl_obj_save, mass, rho, p_atmo, gravity, max_trim, max_heel, cog_ini_ref, cor_stab, z_water_fixed, alfa_fixed, heel)
                #CONTROLLARE hydrostatic_solver("weight_lcg_heel", stl_obj, mass, rho, p_atmo, gravity, max_trim, max_heel, cog_ini_ref, cor_stab, z_water_fixed, alfa_fixed, -heel)

            stl_obj_rotate_x  = rotate_stl_x(stl_obj_save, heel_opt, cog_ini_ref)
            stl_obj_rotate_xy = rotate_stl_y(stl_obj_rotate_x, alfa_opt, cog_ini_ref)
            fz, my, mx = compute_forces_moment(stl_obj_rotate_xy, z_water_opt, cog_ini_ref, rho, p_atmo, gravity)
            # righting arm must be positive if positive mx. however, mx is opposite of real mx
            mx = -mx
            by = mx/fz
            # my should be zero and fz=mass*gravity
            fs.write(f"{myround(heel)} {myround(z_water_opt)} {myround(alfa_opt)} {myround(fz)} {myround(my)} {-myround(mx)} {myround(by)}\n")
        fs.close()
        print_log("...completed")
