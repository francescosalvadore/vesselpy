# trace generated using paraview version 5.5.0

import sys
import os
PARAVIEW_FOLDER = os.environ["PARAVIEW_FOLDER"] 
sys.path.append(PARAVIEW_FOLDER)

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def extract_immersed_hull(filename, filename_out, z_water, include_wl=True):
    # create a new 'STL Reader'
    stl_handle = STLReader(FileNames=[filename])
    
    if include_wl:
        methods = ["clipClosedSurface", "appendSlice"]
        # second method is not able to produce a water-tight immersed hull, and, as a consequence,
        # calculation of center of bouyancy is wrong. kept here, for reference. But only
        # method "0" (clipClosedSurface) must be used.
        method = methods[0]
        if method == methods[0]:
            # create a new 'Clip Closed Surface'
            clipClosedSurface1 = ClipClosedSurface(Input=stl_handle)
            clipClosedSurface1.ClippingPlane = 'Plane'
            clipClosedSurface1.ClippingPlane.Origin = [0.9497355222702026, 0.0, z_water]
            # for some reason this filter normal behaves in the opposite way of the other clip (-1 instead of 1 required)
            clipClosedSurface1.ClippingPlane.Normal = [0.0, 0.0, -1.0]
            # create a new 'Triangulate'
            triangulate1 = Triangulate(Input=clipClosedSurface1)
        elif method == methods[1]:
            clip1 = Clip(Input=stl_handle)
            clip1.ClipType = 'Plane'
            clip1.ClipType.Normal = [0.0, 0.0, 1.0]
            clip1.ClipType.Origin = [0.9497355222702026, 0.0, z_water]
            extractSurface1 = ExtractSurface(Input=clip1)
            triangulate_clip = Triangulate(Input=extractSurface1)

            slice1 = Slice(Input=stl_handle)
            slice1.SliceType = 'Plane'
            slice1.SliceOffsetValues = [0.0]
            slice1.SliceType.Origin = [0.9485052177042235, 0.0, z_water]
            slice1.SliceType.Normal = [0.0, 0.0, 1.0]
            delaunay2D1 = Delaunay2D(Input=slice1) # Create Delaunay
            triangulate_slice = Triangulate(Input=delaunay2D1)

            triangulate1 = AppendGeometry(Input=[triangulate_clip, triangulate_slice])
    else:
        clip1 = Clip(Input=stl_handle)
        clip1.ClipType = 'Plane'
        clip1.ClipType.Normal = [0.0, 0.0, 1.0]
        clip1.ClipType.Origin = [0.9497355222702026, 0.0, z_water]
        extractSurface1 = ExtractSurface(Input=clip1)
        # create a new 'Triangulate'
        triangulate1 = Triangulate(Input=extractSurface1)
    
    # save data
    #new_filename = filename.split(".stl")[0]+"_immersed.stl"
    #new_filename = "hull_immersed.stl"
    new_filename = filename_out
    SaveData(new_filename, proxy=triangulate1, FileType='Ascii')

def extract_wl(filename, filename_out, z_water):
    # create a new 'STL Reader'
    stl_handle = STLReader(FileNames=[filename])

    # Create slice of waterline
    slice1 = Slice(Input=stl_handle)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.SliceType.Origin = [0.9485052177042235, 0.0, z_water]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # Create Delaunay
    delaunay2D1 = Delaunay2D(Input=slice1)

    # Triangulate
    triangulate1 = Triangulate(Input=delaunay2D1)

    new_filename = filename_out

    improve_quality = True
    if improve_quality:
        decimate1 = Decimate(registrationName='Decimate1', Input=triangulate1)
        subdivide1 = Subdivide(registrationName='Subdivide1', Input=decimate1)
        subdivide1.NumberofSubdivisions = 3
        SaveData(new_filename, proxy=subdivide1, FileType='Ascii')
    else:
        SaveData(new_filename, proxy=triangulate1, FileType='Ascii')

if __name__ == "__main__":
    z_water = 0.0912871311
    z_water = 0.08905071944808869
    #extract_immersed_hull(filepaths[0], z_water)
    extract_immersed_hull("hull_rotated.stl", "hull_immersed.stl", z_water)
