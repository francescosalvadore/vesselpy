import math

def read_stl(filepath):
    coordx=[] ; coordy=[] ; coordz=[]
    triangles = [] ; normals = []
    with open(filepath,'r') as f:
        for el in f:
            if len(el.split()) == 4:
                cx = float(el.split()[1]) ; cy = float(el.split()[2]) ; cz = float(el.split()[3])
                coordx.append(cx) ; coordy.append(cy) ; coordz.append(cz)
            if len(el.split()) == 5:
                normals.append([float(el.split()[2]), float(el.split()[3]), float(el.split()[4])])
    #print("nx,ny,nz: ",len(coordx),len(coordy),len(coordz))
    #print("coordx: ",coordx[-1])

    i_tri = -1
    i_co = -1
    for cx, cy, cz in zip(coordx, coordy, coordz):
        i_co += 1
        if i_co%3 == 0:
            i_tri += 1
            triangles.append([])
        triangles[i_tri].append([cx, cy, cz])

    stl_obj = dict(triangles=triangles, normals=normals)

    create_areas_centers_bb(stl_obj)

    return stl_obj

def save_stl(stl_obj, filepath):
    f = open(filepath, "w")
    f.write("solid rotated_hull\n")
    for triangle, center, normal in zip(stl_obj["triangles"], stl_obj["centers"], stl_obj["normals"]):
        t1 = triangle[0]
        t2 = triangle[1]
        t3 = triangle[2]
        f.write(f" facet normal {normal[0]} {normal[1]} {normal[2]}\n")
        f.write(f"  outer loop\n")
        f.write(f"   vertex {t1[0]} {t1[1]} {t1[2]}\n")
        f.write(f"   vertex {t2[0]} {t2[1]} {t2[2]}\n")
        f.write(f"   vertex {t3[0]} {t3[1]} {t3[2]}\n")
        f.write(f"  endloop\n")
        f.write(f" endfacet\n")
    f.write("endsolid\n")
    f.close()

def refine_stl(stl_obj, refine_ratio=10.):
    bb = stl_obj["boundingbox"]
    min_bb = min(bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2])
    max_edge = min_bb/refine_ratio
    print("Refining with max edge: ",max_edge)
    #n_refine_max = 1 
    n_refine_max = 20
    ref_needed = True
    for i in range(n_refine_max):
        if not ref_needed:
            break
        ref_needed = False
        triangles = stl_obj["triangles"]
        print(" --- refine iteration: ",i, " - n_triangles :",len(triangles)) #, end = "")
        normals = stl_obj["normals"]
        triangles_refine = []
        normals_refine = []
        for t,n in zip(triangles, normals):
            p1x = t[0][0] ; p1y = t[0][1] ; p1z = t[0][2]
            p2x = t[1][0] ; p2y = t[1][1] ; p2z = t[1][2]
            p3x = t[2][0] ; p3y = t[2][1] ; p3z = t[2][2]
            l1q = (p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2
            l2q = (p1x-p3x)**2+(p1y-p3y)**2+(p1z-p3z)**2
            l3q = (p2x-p3x)**2+(p2y-p3y)**2+(p2z-p3z)**2
            l1 = math.sqrt(l1q)
            l2 = math.sqrt(l2q)
            l3 = math.sqrt(l3q)
            #RIMETTEREif max(l1,l2,l3) > max_edge:
            tol = 1e-8
            if max(l1,l2,l3) > max_edge and l1q > tol and l2q > tol and l3q > tol:
                # decompose a triangle into three triangles using the barycenter point.
                # beware: order of vertexes is important to keep correct the normal orientation
                bar_x = (p1x+p2x+p3x)/3.
                bar_y = (p1y+p2y+p3y)/3.
                bar_z = (p1z+p2z+p3z)/3.
                t1 = [[p1x,p1y,p1z],[p2x,p2y,p2z],[bar_x,bar_y,bar_z]]
                t2 = [[p1x,p1y,p1z],[bar_x,bar_y,bar_z],[p3x,p3y,p3z]]
                t3 = [[p2x,p2y,p2z],[p3x,p3y,p3z],[bar_x,bar_y,bar_z]]
                triangles_refine.extend([t1,t2,t3])
                normals_refine.extend([n,n,n])
                ref_needed = True
            else:
                triangles_refine.append(t)
                normals_refine.append(n)

        stl_obj["triangles"]=triangles_refine
        stl_obj["normals"]=normals_refine

    create_areas_centers_bb(stl_obj)

def create_areas_centers_bb(stl_obj):
    triangles = stl_obj["triangles"]
    normals = stl_obj["normals"]
    centers = []; areas = []
    minx = 2**20;  miny = 2**20;  minz = 2**20
    maxx = -2**20; maxy = -2**20; maxz = -2**20
    for t,n in zip(triangles, normals):
        #print("t,n: ",t,n)
        p1x = t[0][0] ; p1y = t[0][1] ; p1z = t[0][2]
        p2x = t[1][0] ; p2y = t[1][1] ; p2z = t[1][2]
        p3x = t[2][0] ; p3y = t[2][1] ; p3z = t[2][2]
        minx = min(minx, p1x) ; minx = min(minx, p2x) ; minx = min(minx, p3x)
        miny = min(miny, p1y) ; miny = min(miny, p2y) ; miny = min(miny, p3y)
        minz = min(minz, p1z) ; minz = min(minz, p2z) ; minz = min(minz, p3z)
        maxx = max(maxx, p1x) ; maxx = max(maxx, p2x) ; maxx = max(maxx, p3x)
        maxy = max(maxy, p1y) ; maxy = max(maxy, p2y) ; maxy = max(maxy, p3y)
        maxz = max(maxz, p1z) ; maxz = max(maxz, p2z) ; maxz = max(maxz, p3z)
        l1 = math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2)
        l2 = math.sqrt((p1x-p3x)**2+(p1y-p3y)**2+(p1z-p3z)**2)
        l3 = math.sqrt((p2x-p3x)**2+(p2y-p3y)**2+(p2z-p3z)**2)
        s = (l1+l2+l3)/2
        try:
            area = math.sqrt(s*(s-l1)*(s-l2)*(s-l3))
        except:
            print("Warning! Degenere face Erone negative component found: =",s,s-l1,s-l2,s-l3)
            area = 0.
        cx = t[0][0]+t[1][0]+t[2][0]
        cy = t[0][1]+t[1][1]+t[2][1]
        cz = t[0][2]+t[1][2]+t[2][2]
        c = [cx/3., cy/3., cz/3.]
        centers.append(c)
        areas.append(area)
        #print("t: ",t, c, n, area)
        #break

    total_area = sum(areas)
    #print("Total area: ",total_area)
    #print("Area over bb_yz factor: ",total_area/((maxx-minx)*(maxy-miny)))

    stl_obj["centers"]     = centers
    stl_obj["areas"]       = areas
    stl_obj["boundingbox"] = [minx, miny, minz, maxx, maxy, maxz]
    #print("Bounding box: ")
    #print("minx, miny, minz: ",minx, miny, minz)
    #print("maxx, maxy, maxz: ",maxx, maxy, maxz)

    return

def rotate(origin, point, cos_angle, sin_angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle must be given as cos_angle and sin_angle.
    """
    ox, oy = origin
    px, py = point
    qx = ox + cos_angle * (px - ox) - sin_angle * (py - oy)
    qy = oy + sin_angle * (px - ox) + cos_angle * (py - oy)
    return qx, qy

def rotate_stl_y(stl_obj, alfa, cog):
    triangles = stl_obj["triangles"]
    normals = stl_obj["normals"]
    triangles_rotate = []
    normals_rotate = []
    angle = math.pi/180.*alfa
    cos_alfa = math.cos(angle)
    sin_alfa = math.sin(angle)
    for t,n in zip(triangles, normals):
        #print(t,n)
        triangle = []
        tx = t[0][0] ; ty = t[0][1] ; tz = t[0][2]
        r_tx, r_tz = rotate((cog[0],cog[2]), (tx,tz), cos_alfa, sin_alfa)
        triangle.append([r_tx, ty, r_tz])
        tx = t[1][0] ; ty = t[1][1] ; tz = t[1][2]
        r_tx, r_tz = rotate((cog[0],cog[2]), (tx,tz), cos_alfa, sin_alfa)
        triangle.append([r_tx, ty, r_tz])
        tx = t[2][0] ; ty = t[2][1] ; tz = t[2][2]
        r_tx, r_tz = rotate((cog[0],cog[2]), (tx,tz), cos_alfa, sin_alfa)
        triangle.append([r_tx, ty, r_tz])
        triangles_rotate.append(triangle)
        r_nx, r_nz = rotate((0.,0.), (n[0],n[2]), cos_alfa, sin_alfa)
        normals_rotate.append((r_nx, n[1], r_nz))
    stl_obj_rotate = dict(triangles=triangles_rotate, normals=normals_rotate)
    create_areas_centers_bb(stl_obj_rotate)
    return stl_obj_rotate

def rotate_stl_x(stl_obj, alfa, cog):
    triangles = stl_obj["triangles"]
    normals = stl_obj["normals"]
    triangles_rotate = []
    normals_rotate = []
    angle = math.pi/180.*alfa
    cos_alfa = math.cos(angle)
    sin_alfa = math.sin(angle)
    for t,n in zip(triangles, normals):
        #print(t,n)
        triangle = []
        tx = t[0][0] ; ty = t[0][1] ; tz = t[0][2]
        r_ty, r_tz = rotate((cog[1],cog[2]), (ty,tz), cos_alfa, sin_alfa)
        triangle.append([tx, r_ty, r_tz])
        tx = t[1][0] ; ty = t[1][1] ; tz = t[1][2]
        r_ty, r_tz = rotate((cog[1],cog[2]), (ty,tz), cos_alfa, sin_alfa)
        triangle.append([tx, r_ty, r_tz])
        tx = t[2][0] ; ty = t[2][1] ; tz = t[2][2]
        r_ty, r_tz = rotate((cog[1],cog[2]), (ty,tz), cos_alfa, sin_alfa)
        triangle.append([tx, r_ty, r_tz])
        triangles_rotate.append(triangle)
        r_ny, r_nz = rotate((0.,0.), (n[1],n[2]), cos_alfa, sin_alfa)
        normals_rotate.append((n[0], r_ny, r_nz))
    stl_obj_rotate = dict(triangles=triangles_rotate, normals=normals_rotate)
    create_areas_centers_bb(stl_obj_rotate)
    return stl_obj_rotate
