#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SToGS)
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all gdml)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#
# Add specific modules
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE("${CMAKE_MODULE_PATH}")
#----------------------------------------------------------------------------
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


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/csrc)
include_directories(${CMAKE_BINARY_DIR})
#
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
set(SToGS_SRC
# Detector related
	csrc/SToGS_MaterialConsultant.cc
	csrc/SToGS_DetectorFactory.cc
    csrc/SToGS_RoomsDF.cc
    csrc/SToGS_ScintillatorDF.cc
    csrc/SToGS_MyStoreDF.cc
    csrc/SToGS_LoadFromDetectorFactory.cc
    csrc/SToGS_TwoShellsDetectorConstruction.cc
    csrc/SToGS_GenericsDF.cc
    csrc/SToGS_ArraysDF.cc
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.cc
    csrc/SToGS_ModularPhysicsList.cc
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.cc
# User's Action related
	csrc/SToGS_G4_CopClusterHit.cc
	csrc/SToGS_G4_CopClusterSD.cc
	csrc/SToGS_G4_TrackerHit.cc
	csrc/SToGS_G4_TrackerSD.cc
	csrc/SToGS_UserActionInitialization.cc
    csrc/SToGS_PrintOut.cc
    csrc/SToGS_Ascii.cc
	csrc/SToGS_UserActionManager.cc
)
set(SToGS_HEADER
# Detector related
	csrc/SToGS_MaterialConsultant.hh
	csrc/SToGS_DetectorFactory.hh
    csrc/SToGS_RoomsDF.hh
    csrc/SToGS_ScintillatorDF.hh
    csrc/SToGS_MyStoreDF.hh
    csrc/SToGS_LoadFromDetectorFactory.hh
    csrc/SToGS_TwoShellsDetectorConstruction.hh
    csrc/SToGS_GenericsDF.hh
    csrc/SToGS_ArraysDF.hh
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.hh
    csrc/SToGS_ModularPhysicsList.hh
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.hh
# User's Action related
	csrc/SToGS_G4_CopClusterHit.hh
	csrc/SToGS_G4_CopClusterSD.hh
	csrc/SToGS_G4_TrackerHit.hh
	csrc/SToGS_G4_TrackerSD.hh
	csrc/SToGS_UserActionManager.hh
    csrc/SToGS_PrintOut.hh
    csrc/SToGS_Ascii.hh
	csrc/SToGS_UserActionInitialization.hh
)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# ADD My Plugins if any
#
if ( NOT MYPATH )
    set(MYPATH SToGS_MY)
endif()
if (MYDET)
 if ( EXISTS "${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.hh"
       AND EXISTS "${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.cc" )
	#
	ADD_DEFINITIONS(-DHAS_MYDET)
	#
	set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.cc)
	set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.hh)
    #
    include_directories(${PROJECT_SOURCE_DIR}/${MYPATH})
	#
	MESSAGE(" --> ADD ${MYDET}DetectorConstruction ")
 else()
    MESSAGE(" --> Cannot find " ${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.hh)
    MESSAGE(" --> Or " ${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction.cc)
  endif()
    MESSAGE(" --> Compile with ${MYDET} " ${PROJECT_SOURCE_DIR}/${MYPATH}/${MYDET}DetectorConstruction)
endif()

# ADD MyXXXX
#
# add MyDetectorConstruction if found
if ( EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.hh" AND EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.cc" )
	#
	ADD_DEFINITIONS(-DHAS_MYDET)
	#
	set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.cc) 
	set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/MyDetectorConstruction.hh) 
	#
	MESSAGE(" --> ADD MyDetectorConstruction ")
endif()
if ( EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyPrimaryGeneratorAction.hh" AND EXISTS "${PROJECT_SOURCE_DIR}/csrc/MyPrimaryGeneratorAction.cc")
	#
	ADD_DEFINITIONS(-DHAS_MYPRIMARYGENERATOR)
	#
	set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/MyPrimaryGeneratorAction.cc)
	set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/MyPrimaryGeneratorAction.hh)
	#
	MESSAGE(" --> ADD MyPrimaryGeneratorAction ")
endif()
#
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/csrc/SToGSConfig.h.cmake"
  "${CMAKE_BINARY_DIR}/SToGSConfig.hh"
)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(SToGS_Source SToGS_Source.cc ${SToGS_SRC} ${SToGS_HEADER})
add_executable(SToGS_DetectorBuilder SToGS_DetectorBuilder.cc ${SToGS_SRC} ${SToGS_HEADER})

#
target_link_libraries(SToGS_Source ${Geant4_LIBRARIES} ${EXTRA_LIB})
target_link_libraries(SToGS_DetectorBuilder ${Geant4_LIBRARIES} ${EXTRA_LIB})


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SToGS_Source DESTINATION bin)
install(TARGETS SToGS_DetectorBuilder DESTINATION bin)


