import math

def sort_by_z(triangle):
    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]
    z_values = [v1[2],v2[2],v3[2]]
    i1 = z_values.index(min(z_values))
    z_values[i1] = 2**100
    i2 = z_values.index(min(z_values))
    z_values[i2] = 2**100
    i3 = z_values.index(min(z_values))
    triangle = [triangle[i1], triangle[i2], triangle[i3]]
    return triangle

def segment_z_plane_intersection(segment, zp, tol=1.0e-10):
    v1   = segment[0]
    v2   = segment[1]
    dz_1 = v2[2] - zp
    dz_2 = zp    - v1[2]
    if dz_1 < 0. or dz_2 < 0.:
        print("Error! Segment points must be sorted according to ascending z")
        raise
    if dz_1 + dz_2 < tol:
        scal = 0.
    else:
        scal = dz_1/(dz_1+dz_2)
    x_i = v1[0] + (v2[0]-v1[0])*scal
    y_i = v1[1] + (v2[1]-v1[1])*scal
    z_i = zp
    return [x_i, y_i, z_i]

def area_center_triangle(t):
    try:
        p1x = t[0][0] ; p1y = t[0][1] ; p1z = t[0][2]
    except:
        print("t: ",str(t))
    p2x = t[1][0] ; p2y = t[1][1] ; p2z = t[1][2]
    p3x = t[2][0] ; p3y = t[2][1] ; p3z = t[2][2]
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
    center = [cx/3., cy/3., cz/3.]
    return area, center

def triangle_z_plane_intersection(triangle, zp, tol=1.0e-10):
    triangle = sort_by_z(triangle)
    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]
    if   v1[2] >= zp:
        #print("all above")
        triangles_up   = [triangle]
        triangles_down = []
    elif v3[2] <= zp:
        #print("all below")
        triangles_up   = []
        triangles_down = [triangle]
    elif v2[2] >= zp:
        #print("v2 and v3 above")
        P_13 = segment_z_plane_intersection([v1,v3], zp, tol)
        P_12 = segment_z_plane_intersection([v1,v2], zp, tol)
        triangles_down = [[v1, P_13, P_12]]
        triangles_up   = [[v2, P_13, P_12], [v2, v3,   P_13]]
    elif v2[2] <= zp:
        #print("only v3 above")
        P_13 = segment_z_plane_intersection([v1,v3], zp, tol)
        P_23 = segment_z_plane_intersection([v2,v3], zp, tol)
        triangles_up   = [[v3, P_13, P_23]]
        triangles_down = [[v1, P_13, P_23], [v1, v2,   P_23]]
    else:
        print("Error! Cannot be here")
        raise
    return triangles_up, triangles_down

if __name__ == "__main__":
    p1 = [1.,3.,2.]
    p2 = [1.,5.,4.]
    p3 = [3.,3.,4]
    zp = 3
    triangle = [p1,p2,p3]
    print(f"triangle: {str(triangle)}")
    print(f"plane   : z={zp}")

    triangles_up, triangles_down = triangle_z_plane_intersection(triangle, zp)
    print(f"triangle_down = {triangles_down}")
    print(f"triangle_up   = {triangles_up}")

