#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Paris)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#
# Add specific modules
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

MESSAGE("${CMAKE_MODULE_PATH}")


#----------------------------------------------------------------------------
# Try and Find Root. If found, set HAS_ROOT
#
find_package(ROOT)
if (ROOT_FOUND)
	MESSAGE("ROOT has been found: ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR} ==> Add TTree facilities")
	if (EXISTS ${CMAKE_SOURCE_DIR}/analysis/libParisEvents.so)
		#
		ADD_DEFINITIONS(-DHAS_ROOT)
		# add ROOT includes and libs + ParisEvents (includes and lib which is compiled)
		include_directories(${CMAKE_SOURCE_DIR}/analysis)
		include_directories(${ROOT_INCLUDE_DIR})
		#
		set (EXTRA_LIB ${ROOT_LIBRARIES} -L${CMAKE_SOURCE_DIR}/analysis -lParisEvents)
		MESSAGE(" --> Cool, libParisEvents.so found in analysis !")	
		#
	else()
		MESSAGE(" --> libParisEvents.so NOT FOUND in analysis, you may needed to compile it first to fully exploit ROOT facilities")	
	endif()
else()
	MESSAGE(" --> ROOT has NOT been found")
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/csrc)

if(Geant4_gdml_FOUND)
	ADD_DEFINITIONS(-DHAS_GDML)
	MESSAGE(" --> GDML is ON")	
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
 
set(SToGS_SRC
# Detector related
	csrc/SToGS_MaterialConsultant.cc
	csrc/SToGS_DetectorFactory.cc
    csrc/SToGS_RoomsDF.cc
    csrc/SToGS_ArraysDF.cc
    csrc/SToGS_ScintillatorDF.cc
    csrc/SToGS_LoadFromDetectorFactory.cc
    csrc/SToGS_ShellDetectorConstruction.cc
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.cc
    csrc/SToGS_ModularPhysicsList.cc
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.cc
# User's Action related
	csrc/SToGS_G4_CaloHit.cc
	csrc/SToGS_G4_CaloSD.cc
	csrc/SToGS_G4_SingleHit.cc
	csrc/SToGS_G4_TrackerSD.cc
	csrc/SToGS_UserActionManager.cc
#
#	csrc/SToGS_3by3CubicConstruction.cc
#	csrc/SToGS_Ascii.cc
#	csrc/SToGS_BaseAscii.cc
#	csrc/SToGS_BaseROOT.cc
#	csrc/SToGS_BasicEventAction.cc
#	csrc/SToGS_BasicPrimaryGeneratorAction.cc
#	csrc/SToGS_BasicRunAction.cc
#	csrc/SToGS_BasicSteppingAction.cc
#	csrc/SToGS_CaloHit.cc
#	csrc/SToGS_CaloSD.cc
#	csrc/SToGS_GDRPrimaryGeneratorAction.cc
#	csrc/SToGS_GeneralPhysics.cc
#	csrc/SToGS_LoadDetectorConstruction.cc
#	csrc/SToGS_LowEnergyEMPhysicsList.cc
#	csrc/SToGS_MaterialConsultant.cc
#	csrc/SToGS_OutputManager.cc
#	csrc/SToGS_PenelopeEMPhysicsList.cc
#	csrc/SToGS_PhysicsList.cc
#	csrc/SToGS_PrintOut.cc
#	csrc/SToGS_ROOT.cc
#	csrc/SToGS_SDCRotated.cc
#	csrc/SToGS_SDCRotated_CsI.cc
#	csrc/SToGS_SegmentedDetectorConstruction.cc
#	csrc/SToGS_SegmentedShellDetectorConstruction.cc
#	csrc/SToGS_SemiSphe2LayersConstruction.cc
#	csrc/SToGS_ShellDetectorConstruction.cc
#	csrc/SToGS_SingleHit.cc
#	csrc/SToGS_StandardEMPhysicsList.cc
#	csrc/SToGS_TrackerSD.cc
#	csrc/SToGS_Tree.cc
#	csrc/SToGS_UniformPrimaryGeneratorAction.cc
#	csrc/SToGS_Hadron0.cc
#	csrc/SToGS_Hadron1.cc
#	csrc/SToGS_OpticalROOT.cc
#	csrc/SToGS_BasicTrackingAction.cc
)

set(SToGS_HEADER
# Detector related
	csrc/SToGS_MaterialConsultant.hh
	csrc/SToGS_DetectorFactory.hh
    csrc/SToGS_RoomsDF.hh
    csrc/SToGS_ArraysDF.hh
    csrc/SToGS_ScintillatorDF.hh
    csrc/SToGS_LoadFromDetectorFactory.hh
    csrc/SToGS_ShellDetectorConstruction.hh
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.hh
    csrc/SToGS_ModularPhysicsList.hh
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.hh
# User's Action related
	csrc/SToGS_G4_CaloHit.hh
	csrc/SToGS_G4_CaloSD.hh
	csrc/SToGS_G4_SingleHit.hh
	csrc/SToGS_G4_TrackerSD.hh
	csrc/SToGS_UserActionManager.hh
#
#	csrc/SToGS_3by3CubicConstruction.hh
#	csrc/SToGS_Ascii.hh
#	csrc/SToGS_BaseAscii.hh
#	csrc/SToGS_BaseROOT.hh
#	csrc/SToGS_BasicEventAction.hh
#	csrc/SToGS_BasicPrimaryGeneratorAction.hh
#	csrc/SToGS_BasicRunAction.hh
#	csrc/SToGS_BasicSteppingAction.hh
#	csrc/SToGS_CaloHit.hh
#	csrc/SToGS_CaloSD.hh
#	csrc/SToGS_GDRPrimaryGeneratorAction.hh
#	csrc/SToGS_GeneralPhysics.hh
#	csrc/SToGS_IonPhysList.hh
#	csrc/SToGS_IonPrimaryGeneratorAction.hh
#	csrc/SToGS_LoadDetectorConstruction.hh
#	csrc/SToGS_LowEnergyEMPhysicsList.hh
#	csrc/SToGS_MaterialConsultant.hh
#	csrc/SToGS_OutputManager.hh
#	csrc/SToGS_PenelopeEMPhysicsList.hh
#	csrc/SToGS_PhysicsList.hh
#	csrc/SToGS_PrintOut.hh
#	csrc/SToGS_ROOT.hh
#	csrc/SToGS_SDCRotated.hh
#	csrc/SToGS_SDCRotated_CsI.hh
#	csrc/SToGS_SegmentedDetectorConstruction.hh
#	csrc/SToGS_SegmentedShellDetectorConstruction.hh
#	csrc/SToGS_SemiSphe2LayersConstruction.hh
#	csrc/SToGS_ShellDetectorConstruction.hh
#	csrc/SToGS_SingleHit.hh
#	csrc/SToGS_StandardEMPhysicsList.hh
#	csrc/SToGS_TrackerSD.hh
#	csrc/SToGS_Tree.hh
#	csrc/SToGS_UniformPrimaryGeneratorAction.hh
#	csrc/SToGS_Hadron0.hh
#	csrc/SToGS_Hadron1.hh
#	csrc/SToGS_Hadron1.hh
#	csrc/SToGS_OpticalROOT.hh
#	csrc/SToGS_BasicTrackingAction.hh
)

# add MyDetectorConstruction if found
if ( EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.hh" AND EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.cc" )
	#
	ADD_DEFINITIONS(-DHAS_MYDETECTOR)
	#
	set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.cc) 
	set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.hh) 
	#
	MESSAGE(" --> ADD MyDetectorConstruction, # 100 in global.paris")
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#add_executable(ParisSource ParisSource.cc ${SToGS_SRC} ${SToGS_HEADER})
add_executable(SToGS_DetectorBuilder SToGS_DetectorBuilder.cc ${SToGS_SRC} ${SToGS_HEADER})

#
#target_link_libraries(ParisSource ${Geant4_LIBRARIES} ${EXTRA_LIB})
target_link_libraries(SToGS_DetectorBuilder ${Geant4_LIBRARIES} ${EXTRA_LIB})


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS ParisSource DESTINATION bin)
install(TARGETS SToGS_DetectorBuilder DESTINATION bin)


