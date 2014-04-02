# strings are recognized when not commented by '#' and their values is assumed to be alone on the next line with no comments in between !
#
# cube : first scintillator ; cuboid : second scintillator
# OR
# cuboid : unique scintillator
#
# 'whichGeometry'
# IS MANDATORY AND NEEDS TO BE FIRST - both with the highest z face equipped with a photocathode :
#	0 : lonely cuboid (only 'cuboid' parameters will be read, except for the cube_side)
#	1 : cube/cuboid assembly
# 'whichPhotocathodeMat' : 0 for aluminium (working), 1 for SbKCs (density not correct yet)
# 'cube_side' : for a PW, the size of the cube ; for a lonely cuboid, the size of x and y sides (in inch)
# 'cuboid_side' : the size of the cuboid z (second scintillator for a PW or z for the lonely cuboid) (in inch)
# 'cube_position' : position of the cube center in cm
# 'cuboid_position' : position of the cuboid center in cm
# 'scint_cube_reflectivity' : reflectivity between the world and the cube ; 2 doubles are expected (it is an array)
# 'scint_cuboid_reflectivity' : reflectivity between the world and the cuboid ; 2 doubles are expected (it is an array)
# 'cube_material' : LaBr3, NaI or CsI
# 'cuboid_material' : LaBr3, NaI or CsI
# 'whichOpticalSurface' : (setModel, setType, setFinish) tests ; see customDetector.cc for more details
#	0 : Specular reflection
#	4 : Teflon (don't work for NaI and CsI)
#	5 : Lambertian reflection
# 'cube_cuboid_rindex' : array (2 doubles expected) for the refractive index between the cube and the cuboid (NO BORDER DEFINED FOR NOW = not read parameter)
# 'cube_cuboid_reflectivity' : reflectivity between the cube and the cuboid ; 2 doubles are expected (it is an array) (NO BORDER DEFINED FOR NOW = not read parameter)
# 'efficiencies' : changes nothing in current model ; ignored for now
#
################################################################################
############################ NaI PHOSWICH EXAMPLE ##############################
############################     2"X2"X(2"+6")    ##############################
############################ Lambertian reflector ##############################
################################################################################
#
whichGeometry
1
#
whichPhotocathMat
0
#
cube_side
2.0
#
cuboid_side
6.0
#
cube_position
0.0
#
cuboid_position
10.16
#
scint_cube_reflectivity
0.97
0.97
#
scint_cuboid_reflectivity
0.97
0.97
#
#cube_cuboid_reflectivity
#0.0
#0.0
#
cube_material
LaBr3
#
cuboid_material
NaI
#
whichOpticalSurface_cube
5
#
whichOpticalSurface_cuboid
5
#
#whichOpticalSurface_cube_cuboid
#0
#
#scint_cube_efficiency
#0.0
#0.0
#
#scint_cuboid_efficiency
#0.0
#0.0
#
#cube_cuboid_efficiency
#0.0
#0.0
#
#cube_cuboid_rindex
#3.0
#3.0