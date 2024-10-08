from scipy.optimize import root, root_scalar
#from scipy.optimize import fsolve
from triangle_utils import triangle_z_plane_intersection, area_center_triangle
from pv_base import extract_immersed_hull
from of_base import get_features_from_of
from stl_lib import read_stl, refine_stl, rotate_stl_y, rotate_stl_x, save_stl, rotate
from utils import print_log
import math

def compute_forces_moment(stl_obj, z_water, cor, rho, p_atmo, gravity):
    fx, fz, my = 0., 0., 0.
    wet_area = 0.
    mx, fy =  0., 0.
    for triangle, center, normal, area in zip(stl_obj["triangles"], stl_obj["centers"],stl_obj["normals"], stl_obj["areas"]):
        triangles_up, triangles_down = triangle_z_plane_intersection(triangle, z_water)
        for t_up in triangles_up:
            area, center = area_center_triangle(t_up)
            x, y, z = center
            fx_temp = - p_atmo*normal[0]*area
            fx += fx_temp
            fy_temp = - p_atmo*normal[1]*area
            fy += fy_temp
            fz_temp = - p_atmo*normal[2]*area
            fz += fz_temp
            my += fz_temp*(x-cor[0]) - fx_temp*(z-cor[2])
            mx += fz_temp*(y-cor[1]) - fy_temp*(z-cor[2])
        for t_down in triangles_down:
            area, center = area_center_triangle(t_down)
            x, y, z = center
            fx_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[0]*area
            fx += fx_temp
            fy_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[1]*area
            fy += fy_temp
            fz_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[2]*area
            fz += fz_temp
            my += fz_temp*(x-cor[0]) - fx_temp*(z-cor[2])
            mx += fz_temp*(y-cor[1]) - fy_temp*(z-cor[2])
            wet_area += area
    #print_log("Iterating values  : fx, fz, my, wet_area: ",fx, fz, my, wet_area)
    #print_log("Iterating values/2: fy, mx, bx: ",fy, mx, mx/fz)
    #print_log(f"fz={fz}")
    #print_log("...", end="", flush=True)
    return fz, my, mx

def old_compute_forces_moment(stl_obj, z_water, cog, rho, p_atmo, gravity):
    # this version does not split triangles and then may require refinement to achieve good results
    fx, fz, my = 0., 0., 0.
    wet_area = 0.
    mx, fy =  0., 0.
    for triangle, center, normal, area in zip(stl_obj["triangles"], stl_obj["centers"],stl_obj["normals"], stl_obj["areas"]):
        x, y, z = center
        if z < z_water:
            fx_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[0]*area
            fx += fx_temp
            fy_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[1]*area
            fy += fy_temp
            fz_temp = - (p_atmo+rho*gravity*(z_water-z))*normal[2]*area
            fz += fz_temp
            my += fz_temp*(x-cog[0]) - fx_temp*(z-cog[2])
            mx += fz_temp*(y-cog[1]) - fy_temp*(z-cog[2])
            wet_area += area
        else:
            fx_temp = - p_atmo*normal[0]*area
            fx += fx_temp
            fy_temp = - p_atmo*normal[1]*area
            fy += fy_temp
            fz_temp = - p_atmo*normal[2]*area
            fz += fz_temp
            my += fz_temp*(x-cog[0]) - fx_temp*(z-cog[2])
            mx += fz_temp*(y-cog[1]) - fy_temp*(z-cog[2])
    #print_log("Iterating values  : fx, fz, my, wet_area: ",fx, fz, my, wet_area)
    #print_log("Iterating values/2: fy, mx, bx: ",fy, mx, mx/fz)
    print_log("...", end="", flush=True)
    return fz, my, mx

def compute_floatation(stl_wl_obj, z_water_opt):
    lcf = 0.
    tcf = 0.
    area_tot = 0.
    for t, center in zip(stl_wl_obj["triangles"], stl_wl_obj["centers"]):
        cx, cy, _ = center
        p1x = t[0][0] ; p1y = t[0][1] ; p1z = 0.
        p2x = t[1][0] ; p2y = t[1][1] ; p2z = 0.
        p3x = t[2][0] ; p3y = t[2][1] ; p3z = 0.
        l1 = math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2)
        l2 = math.sqrt((p1x-p3x)**2+(p1y-p3y)**2+(p1z-p3z)**2)
        l3 = math.sqrt((p2x-p3x)**2+(p2y-p3y)**2+(p2z-p3z)**2)
        s = (l1+l2+l3)/2
        try:
            area = math.sqrt(s*(s-l1)*(s-l2)*(s-l3))
        except:
            print_log("Warning! Degenere face Erone negative component found: =",s,s-l1,s-l2,s-l3)
            area = 0.
        area_tot += area
        lcf      += cx*area
        tcf      += cy*area
    lcf /= area_tot
    tcf /= area_tot

    mom_in_x = 0.
    mom_in_y = 0.
    for t, center in zip(stl_wl_obj["triangles"], stl_wl_obj["centers"]):
        cx, cy, _ = center
        p1x = t[0][0] ; p1y = t[0][1] ; p1z = 0.
        p2x = t[1][0] ; p2y = t[1][1] ; p2z = 0.
        p3x = t[2][0] ; p3y = t[2][1] ; p3z = 0.
        l1 = math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2)
        l2 = math.sqrt((p1x-p3x)**2+(p1y-p3y)**2+(p1z-p3z)**2)
        l3 = math.sqrt((p2x-p3x)**2+(p2y-p3y)**2+(p2z-p3z)**2)
        s = (l1+l2+l3)/2
        try:
            area = math.sqrt(s*(s-l1)*(s-l2)*(s-l3))
        except:
            print_log("Warning! Degenere face Erone negative component found: =",s,s-l1,s-l2,s-l3)
            area = 0.
        mom_in_x += ((cx-lcf)**2)*area
        mom_in_y += ((cy-tcf)**2)*area

    xmin, ymin, zmin, xmax, ymax, zmax = stl_wl_obj["boundingbox"]
    LWL = xmax - xmin
    BWL = ymax - ymin
    return area_tot, lcf, tcf, LWL, BWL, xmin, xmax, ymin, ymax, mom_in_x, mom_in_y

def compute_variations(stl_obj, z_water_opt, cog, mass, LWL, LOA, rho, p_atmo, gravity):

    z_water_opt_sinked_plus1 = z_water_opt + 0.01
    fz, my, mx = compute_forces_moment(stl_obj, z_water_opt_sinked_plus1, cog, rho, p_atmo, gravity)
    weight_to_sink = fz/gravity-mass
    
    # moment to trim based on different angles is possible. 
    # minus sign is to have positive my. division by gravity to have results in kgf m / cm
    #alfa_1 = - 1.0
    alfa_1 = - math.atan(0.01/LWL)*180/math.pi 
    #alfa_1  = - math.atan(0.01/LOA)*180/math.pi 
    #alfa_1  = - math.atan(0.01/109.027)*180/pi 
    stl_obj_rotate_y = rotate_stl_y(stl_obj, alfa_1, cog)
    fz, my, mx = compute_forces_moment(stl_obj_rotate_y, z_water_opt, cog, rho, p_atmo, gravity)
    moment_to_trim = my/gravity

    # il prossimo e' diviso per il delta z (quindi momento per creare dz=1)
    #alfa_1  = - 1.0
    #delta_z_cm = LWL*math.tan(abs(alfa_1)*pi/180)*100
    ##delta_z_cm = 109.027*math.tan(abs(alfa_1)*pi/180)*100
    #stl_obj_rotate_y = rotate_stl_y(stl_obj, alfa_1, cog)
    #fz, my, mx = compute_forces_moment(stl_obj_rotate_y, z_water_opt, cog, rho, p_atmo, gravity)
    #moment_to_trim = my/(gravity*delta_z_cm)

    return weight_to_sink, moment_to_trim

def equations_3DOF(p, stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_trim, max_heel):
    z_water, alfa, heel = p
    print_log(f"z_water={z_water} alfa={alfa} heel={heel}")
    # must limit alfa, otherwise there are around +/-90 degrees solutions
    if abs(heel) > max_heel:
        return 1000., 1000., 1000.
    if abs(alfa) > max_trim:
        return 1000., 1000., 1000
    # use -heel so that cor_y > 0 corresponds to positive heels
    stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel, cor)
    stl_obj_rotate_xy = rotate_stl_y(stl_obj_rotate_x, alfa, cor)
    fz, my, mx = compute_forces_moment(stl_obj_rotate_xy, z_water, cor, rho, p_atmo, gravity)
    cog_rotate = full_rotate_point(cor, cog, alfa, heel, first_rotation="heel")
    return fz-mass*gravity, my-mass*gravity*(cog_rotate[0]-cor[0]), mx-mass*gravity*(cog_rotate[1]-cor[1])
    #return fz-mass*gravity, my, mx

def equations_2DOF_a(p, stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_trim):
    z_water, alfa = p
    # must limit alfa, otherwise there are around +/-90 degrees solutions
    if abs(alfa) > max_trim:
        return 1000., 1000.
    stl_obj_rotate = rotate_stl_y(stl_obj, alfa, cor)
    fz, my, _ = compute_forces_moment(stl_obj_rotate, z_water, cor, rho, p_atmo, gravity)
    cog_rotate = full_rotate_point(cor, cog, alfa, 0., first_rotation="trim")
    return fz-mass*gravity, my-mass*gravity*(cog_rotate[0]-cor[0])

def equations_2DOF_b_old(p, stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_heel):
    z_water, heel = p
    if abs(heel) > max_heel:
        return 1000., 1000.
    stl_obj_rotate = rotate_stl_x(stl_obj, heel, cor)
    fz, _ , mx = compute_forces_moment(stl_obj_rotate, z_water, cor, rho, p_atmo, gravity)

    #BUGGONE SCAMBIATI HEEL AND 0 cog_rotate = full_rotate_point(cor, cog, heel, 0., first_rotation="heel")
    cog_rotate = full_rotate_point(cor, cog, 0., heel, first_rotation="heel")

    print_log(f"z_water = {z_water} - fz = {fz} - mass*gravity={mass*gravity} - residual = {fz-mass*gravity}/{mx-mass*gravity*(cog_rotate[1]-cor[1])}")
    return fz-mass*gravity, mx-mass*gravity*(cog_rotate[1]-cor[1])

def equations_2DOF_b(p, stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_heel, alfa):
    z_water, heel = p
    if abs(heel) > max_heel:
        return 1000., 1000.
    stl_obj_rotate_x = rotate_stl_x(stl_obj, heel, cor)
    stl_obj_rotate = rotate_stl_y(stl_obj_rotate_x, alfa, cor)
    fz, _ , mx = compute_forces_moment(stl_obj_rotate, z_water, cor, rho, p_atmo, gravity)
    cog_rotate = full_rotate_point(cor, cog, alfa, heel, first_rotation="heel")
    return fz-mass*gravity, mx-mass*gravity*(cog_rotate[1]-cor[1])

def equations_1DOF(p, stl_obj, cor, mass, rho, p_atmo, gravity):
    z_water,  = p
    fz, _ , _ = compute_forces_moment(stl_obj, z_water, cor, rho, p_atmo, gravity)
    print_log(f"z_water = {z_water} - fz = {fz} - mass*gravity={mass*gravity} - residual = {fz-mass*gravity}")
    return fz-mass*gravity, 

def hydrostatic_solver(mode, stl_obj, mass, rho, p_atmo, gravity, max_trim, max_heel, cog, cor, z_water, alfa, heel):
    # Next part should be done later because if a DOF is fixed, rotation is done around provisional cog and bounding box 
    # changes. For small angles not a problem, though. TODO. Old comment, probably not true any more.
    z_water_start = z_water
    alfa_start    = alfa
    heel_start    = heel

    equilibrium_found = False
    if mode == "weight_lcg_tcg":
        #(z_water_opt, alfa_opt, heel_opt), infodict, ierr, msg = fsolve(func=equations_3DOF, x0=(z_water_start, alfa_start, heel_start), args=(stl_obj, cor, mass), full_output=True)
        solution = root(fun=equations_3DOF, x0=(z_water_start, alfa_start, heel_start), args=(stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_trim, max_heel), method="hybr") #"anderson")
        (z_water_opt, alfa_opt, heel_opt) = solution.x
        msg                               = solution.message
        ierr                              = 1 if solution.success else 0
    elif mode == "weight_lcg_heel":
        heel_opt = heel
        stl_obj_rotate_x = rotate_stl_x(stl_obj, heel, cor)
        cog_rotate = full_rotate_point(cor, cog, 0., heel, first_rotation="heel")
        solution = root(fun=equations_2DOF_a, x0=(z_water_start, alfa_start), args=(stl_obj_rotate_x, cog_rotate, cor, mass, rho, p_atmo, gravity, max_trim), method="hybr")
        (z_water_opt, alfa_opt) = solution.x
        msg                     = solution.message
        ierr                    = 1 if solution.success else 0
    elif mode == "weight_trim_tcg":
        alfa_opt = alfa

        # Since the imposed rotation (alfa) is done after heel and heel is unknown, alfa-rotation cannot be done here at the beginning but
        # has to be done during iterations (following right order heel -> alfa)
        #INCORRECT stl_obj_rotate_y = rotate_stl_y(stl_obj, alfa, cor)
        #INCORRECT cog_rotate = full_rotate_point(cor, cog, alfa, 0., first_rotation="heel")
        #INCORRECT solution = root(fun=equations_2DOF_b_old, x0=(z_water_start, heel_start), args=(stl_obj_rotate_y, cog_rotate, cor, mass, rho, p_atmo, gravity, max_heel), method="broyden1")
        solution = root(fun=equations_2DOF_b, x0=(z_water_start, heel_start), args=(stl_obj, cog, cor, mass, rho, p_atmo, gravity, max_heel, alfa), method="hybr")

        (z_water_opt, heel_opt) = solution.x
        msg                     = solution.message
        ierr                    = 1 if solution.success else 0
    elif mode == "weight_trim_heel":
        alfa_opt = alfa
        heel_opt = heel
        stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel, cor)
        stl_obj_rotate_xy = rotate_stl_y(stl_obj_rotate_x, alfa, cor)
        #(z_water_opt,), infodict, ierr, msg = fsolve(func=equations_1DOF, x0=(z_water_start,), args=(stl_obj_rotate_xy, cor, mass), full_output=True)
        solution = root(fun=equations_1DOF, x0=(z_water_start,), args=(stl_obj_rotate_xy, cor, mass, rho, p_atmo, gravity), method="hybr")
        (z_water_opt, ) = solution.x
        msg             = solution.message
        ierr            = 1 if solution.success else 0
    elif mode == "sinkage_trim_heel":
        z_water_opt = z_water
        alfa_opt    = alfa
        heel_opt    = heel
        ierr        = 1

    ret = [z_water_opt, alfa_opt, heel_opt]

    stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel_opt, cor)
    stl_obj_rotate_xy = rotate_stl_y(stl_obj_rotate_x, alfa_opt, cor)

    fz, my, mx = compute_forces_moment(stl_obj_rotate_xy, z_water_opt, cor, rho, p_atmo, gravity)

    mass_recompute  = fz/gravity

    cog_x_recompute = my/fz + cor[0]
    #print_log("corregere qui!!")
    ####return fz-mass*gravity, my-mass*gravity*(cog_rotate[0]-cor[0]), mx-mass*gravity*(cog_rotate[1]-cor[1])
    print_log("cor[0],-my/fz,cog_x_recompute: ",cor[0],-my/fz,cog_x_recompute)

    # the sign of next mx/fz has to be verified!!!
    cog_y_recompute = mx/fz + cor[1]
    #BOHcog_y_recompute = -mx/fz + cor[1]

    # We have just found cog_x/y_recompute in the hydrostatics position, i.e. in the new position/reference of the hull.
    # cog_z_is given in input in the initial reference. Then we have to find cog_z_recompute in the hydrostatics reference
    # such that the corresponding initial-reference cog_z is the given one. We need to solve the impicit problem.
    def equation_cog_z(cog_z_current, cog_z_ini_ref, cog_x, cog_y, cor, alfa_opt, heel_opt):
        cog_ini_ref = full_rotate_point(cor, [cog_x, cog_y, cog_z_current], -alfa_opt, -heel_opt, first_rotation="trim")
        return cog_ini_ref[2]-cog_z_ini_ref

    sol = root_scalar(f=equation_cog_z, x0=cog[2], method="bisect", bracket=[cog[2]-1000,cog[2]+1000], \
                      args=(cog[2], cog_x_recompute, cog_y_recompute, cor, alfa_opt, heel_opt))
    cog_z_recompute = sol.root

    cog_rotate = full_rotate_point(cor, [cog_x_recompute, cog_y_recompute, cog_z_recompute], -alfa_opt, -heel_opt, first_rotation="heel")
    print_log("check cog_rotate: ",cog_rotate)
    ret += [mass_recompute, cog_x_recompute, cog_y_recompute, cog_z_recompute]

    if ierr == 1: 
        equilibrium_found = True
    else:
        equilibrium_found = False
        print_log("Error! Non-linear solver error message: ",msg)
        raise

    return ret+[equilibrium_found]

def get_features_base(stl_obj_save, stl_obj_static, z_water_opt):
    # Overall dimensions
    xmin, ymin, zmin, xmax, ymax, zmax = stl_obj_save["boundingbox"]
    #print_log(f"xmin={xmin} xmax={xmax} ymin={ymin} ymax={ymax} zmin={zmin} zmax={zmax}")
    LOA = xmax - xmin
    BOA = ymax - ymin
    D   = zmax - zmin
    LOA_over_BOA = LOA/BOA
    BOA_over_D   = BOA/D
    #print_log(f"LOA={LOA} - BOA={BOA} - D={D}")
    _, _, zmin_static, _, _, _ = stl_obj_static["boundingbox"]
    T   = z_water_opt - zmin_static
    return LOA, BOA, D, T, xmin, ymin, zmin, xmax, ymax, zmax

def get_form_coefficients(LWL, BWL, T, displaced_volume, AWP, wet_area):
    Cb  = displaced_volume/(LWL*BWL*T)
    Cvp = displaced_volume/(AWP*T)
    Cwp = AWP/(LWL*BWL)
    Cws = wet_area/(displaced_volume*LWL)**0.5
    #Cws = wet_area/(displaced_volume)**(2/3)  alternative form to avoid LWL inaccuracies
    return Cb, Cvp, Cwp, Cws

def get_intact_stability_heel(stl_obj, z_water_opt, cog, cob, cof, rho, p_atmo, gravity, heel_small=0.1):
    # Works on resultant hydrostatics conditions (assumes COG and COB along the same vertical)
    rotation_center = "cog"
    if rotation_center == "ori":
        stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel_small, [0.,0.,0.])
    elif rotation_center == "cog":
        stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel_small, cog)
    elif rotation_center == "cof":
        stl_obj_rotate_x  = rotate_stl_x(stl_obj, heel_small, cof)
    fz, _, mx = compute_forces_moment(stl_obj_rotate_x, z_water_opt, cog, rho, p_atmo, gravity)
    by = mx/fz
    GMt = by/math.sin(math.pi/180*heel_small)
    BMt = GMt + cog[2] - cob[2]
    Mt  = GMt + cog[2] - z_water_opt
    metacenter = [cog[0], cog[1], cog[2]+GMt]
    #print_log(f"heel by={by} Mx={mx} GMt={GMt} metacenter={str(metacenter)}")
    return GMt, BMt, Mt, metacenter

def get_intact_stability_trim(stl_obj, z_water_opt, cog, cob, cof, rho, p_atmo, gravity, trim_small=0.1):
    # Works on resultant hydrostatics conditions (assumes COG and COB along the same vertical)
    rotation_center = "cog"
    if rotation_center == "ori":
        stl_obj_rotate_y  = rotate_stl_y(stl_obj, -trim_small, [0.,0.,0.])
    elif rotation_center == "cog":
        stl_obj_rotate_y  = rotate_stl_y(stl_obj, -trim_small, cog)
    elif rotation_center == "cof":
        stl_obj_rotate_y  = rotate_stl_y(stl_obj, -trim_small, cof)

    fz, my, _ = compute_forces_moment(stl_obj_rotate_y, z_water_opt, cog, rho, p_atmo, gravity)
    bx = my/fz
    GMl = bx/math.sin(math.pi/180*trim_small)
    BMl = GMl + cog[2] - cob[2]
    Ml  = GMl + cog[2] - z_water_opt
    metacenter = [cog[0], cog[1], cog[2]+GMl]
    #print_log(f"trim bx={bx} My={my} GMl={GMl} metacenter={str(metacenter)}")
    return GMl, BMl, Ml, metacenter

def get_metacenter_heel_general(stl_obj, z_water_opt, mass, cog, heel_small=0.1):
    # Get heeling metacenter considering generic algorithm, i.e., intersection of vertical from
    # Center Of Bouyancy considering two conditions: line passing per COG of heel-rotated condition,
    # and line with rotation of small heel angle and rotationvessel reference.
    save_stl(stl_obj, "hull_save.stl")
    extract_immersed_hull("hull_save.stl", "hull_save_immersed.stl", z_water_opt)
    displaced_volume_initial, cob_initial_oldref, _ = get_features_from_of("hull_save_immersed.stl", mass, cog)
    angle = -math.pi/180.*heel_small ; cos_angle = math.cos(angle) ; sin_angle = math.sin(angle)
    cob_initial_temp = rotate([cog[1],cog[2]], [cob_initial_oldref[1],cob_initial_oldref[2]], cos_angle, sin_angle)
    cob_initial = [cob_initial_oldref[0], cob_initial_temp[0], cob_initial_temp[1]]
    stl_obj_save_rotated  = rotate_stl_x(stl_obj, heel_small, cog)
    save_stl(stl_obj_save_rotated, "hull_save_rotated.stl")
    extract_immersed_hull("hull_save_rotated.stl", "hull_save_rotated_immersed.stl", z_water_opt)
    displaced_volume_initial_rotated, cob_initial_rotated, _ = get_features_from_of("hull_save_rotated_immersed.stl", mass, cog)
    #print_log(f"cob_initial={str(cob_initial)} - cob_initial_rotated={str(cob_initial_rotated)}")
    metacenter_x = cob_initial[0]
    metacenter_y = cob_initial_rotated[1]
    metacenter_z = cob_initial[2]+math.tan((90-heel_small)*math.pi/180.)*(cob_initial_rotated[1]-cob_initial[1])
    #print_log(f"heel metacenter={metacenter_x} - metacenter={metacenter_y} - metacenter={metacenter_z}")
    return [metacenter_x, metacenter_y, metacenter_z]

def get_metacenter_trim_general(stl_obj, z_water_opt, mass, cog, trim_small=0.1):
    # Get trimming metacenter considering generic algorithm, i.e., intersection of vertical from
    # Center Of Bouyancy considering two conditions: line passing per COG of heel-rotated condition,
    # and line with rotation of small heel angle and rotationvessel reference.
    save_stl(stl_obj, "hull_save.stl")
    extract_immersed_hull("hull_save.stl", "hull_save_immersed.stl", z_water_opt)
    displaced_volume_initial, cob_initial_oldref, _ = get_features_from_of("hull_save_immersed.stl", mass, cog)
    #print_log(f"TRIMGENERAL cob_initial_oldref={cob_initial_oldref}")
    angle = math.pi/180.*trim_small ; cos_angle = math.cos(angle) ; sin_angle = math.sin(angle)
    cob_initial_temp = rotate([cog[0],cog[2]], [cob_initial_oldref[0],cob_initial_oldref[2]], cos_angle, sin_angle)
    cob_initial = [cob_initial_temp[0], cob_initial_oldref[1], cob_initial_temp[1]]
    stl_obj_save_rotated  = rotate_stl_y(stl_obj, trim_small, cog)
    save_stl(stl_obj_save_rotated, "hull_save_rotated.stl")
    extract_immersed_hull("hull_save_rotated.stl", "hull_save_rotated_immersed.stl", z_water_opt)
    displaced_volume_initial_rotated, cob_initial_rotated, _ = get_features_from_of("hull_save_rotated_immersed.stl", mass, cog)
    #print_log(f"cob_initial={str(cob_initial)} - cob_initial_rotated={str(cob_initial_rotated)}")
    metacenter_x = cob_initial_rotated[0]
    metacenter_y = cob_initial[1]
    metacenter_z = cob_initial[2]-math.tan((90-trim_small)*math.pi/180.)*(cob_initial_rotated[0]-cob_initial[0])
    #print_log(f"trim metacenter={metacenter_x} - metacenter={metacenter_y} - metacenter={metacenter_z}")
    return [metacenter_x, metacenter_y, metacenter_z]

def rotate_point(center, p, angle): #, trim_or_heel):    
    cos_angle  = math.cos(math.pi/180.*angle) ; sin_angle = math.sin(math.pi/180.*angle)
    p_1, p_2 = rotate((center[0],center[1]), (p[0],p[1]), cos_angle, sin_angle)
    return p_1, p_2

def full_rotate_point(center, point, alfa, heel, first_rotation="heel"):
    if   first_rotation == "trim":
        point_rot_x, point_rot_z = rotate_point((center[0],center[2]), (point[0],point[2]),    alfa)
        point_rot_y, point_rot_z = rotate_point((center[1],center[2]), (point[1],point_rot_z), heel)
        return [point_rot_x, point_rot_y, point_rot_z]
    elif first_rotation == "heel":
        point_rot_y, point_rot_z = rotate_point((center[1],center[2]), (point[1],point[2]),    heel)
        point_rot_x, point_rot_z = rotate_point((center[0],center[2]), (point[0],point_rot_z), alfa)
        return [point_rot_x, point_rot_y, point_rot_z]
