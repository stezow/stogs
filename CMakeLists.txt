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
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
MESSAGE("[SToGS] Search for cmake modules in ${CMAKE_MODULE_PATH}")
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/csrc)
include_directories(${CMAKE_BINARY_DIR})
#
set(SToGS_HEADER)
set(SToGS_SRC)
set(SToGS_EXTRA_LIB "")
#
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#
# Add third party libraries (model: -L/path -ltoto ) and includes
#
#
if (THIRD_PARTY_LIBS)
    set (SToGS_EXTRA_LIB ${SToGS_EXTRA_LIB} ${THIRD_PARTY_LIBS})
    MESSAGE("[SToGS] Add third party paths to libraries: " ${THIRD_PARTY_LIBS})
endif()
if (THIRD_PARTY_INCLUDE)
    include_directories(${THIRD_PARTY_INCLUDE})
    MESSAGE("[SToGS] Add third party paths to includes: " ${THIRD_PARTY_INCLUDE})
endif()
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
set(SToGS_SRC ${SToGS_SRC}
# Detector related
	csrc/SToGS_MaterialConsultant.cc
	csrc/SToGS_DetectorFactory.cc
    csrc/SToGS_RoomsDF.cc
    csrc/SToGS_ScintillatorDF.cc
    csrc/SToGS_SemiConductorGeDF.cc
    csrc/SToGS_MyStoreDF.cc
    csrc/SToGS_LoadFromDetectorFactory.cc
    csrc/SToGS_TwoShellsDetectorConstruction.cc
        csrc/SToGS_AGATA.cc
    csrc/SToGS_GenericsDF.cc
    csrc/SToGS_ArraysDF.cc
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.cc
    csrc/SToGS_ModularPhysicsList.cc
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.cc
# User's Action related
	csrc/SToGS_G4_TrackInformation.cc
	csrc/SToGS_G4_CopClusterHit.cc
	csrc/SToGS_G4_CopClusterSD.cc
	csrc/SToGS_G4_TrackerHit.cc
	csrc/SToGS_G4_TrackerSD.cc
	csrc/SToGS_UserActionInitialization.cc
    csrc/SToGS_PrintOut.cc
    csrc/SToGS_Ascii.cc
	csrc/SToGS_UserActionManager.cc
)
set(SToGS_HEADER ${SToGS_HEADER}
# Detector related
	csrc/SToGS_MaterialConsultant.hh
	csrc/SToGS_DetectorFactory.hh
    csrc/SToGS_RoomsDF.hh
    csrc/SToGS_ScintillatorDF.hh
    csrc/SToGS_SemiConductorGeDF.hh
    csrc/SToGS_MyStoreDF.hh
    csrc/SToGS_LoadFromDetectorFactory.hh
    csrc/SToGS_TwoShellsDetectorConstruction.hh
    csrc/SToGS_GenericsDF.hh
    csrc/SToGS_ArraysDF.hh
        csrc/SToGS_AGATA.hh
# Physics Lists related
    csrc/SToGS_G4_GeneralPhysics.hh
    csrc/SToGS_ModularPhysicsList.hh
# Generator related
    csrc/SToGS_G4_GPSPrimaryGeneratorAction.hh
# User's Action related
	csrc/SToGS_G4_TrackInformation.hh
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
# Try and Find Root. If found, set HAS_ROOT in SToGSConfig and compile Core ROOT based modules
#
find_package(ROOT)
if (ROOT_FOUND)
	MESSAGE("[SToGS] ROOT has been found: ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR}")
    #
    ADD_DEFINITIONS(-DHAS_ROOT)
    include_directories(${ROOT_INCLUDE_DIR})
    set (SToGS_EXTRA_LIB ${SToGS_EXTRA_LIB} ${ROOT_LIBRARIES})
    #
    set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOT.hh)
    set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOT.cc)
 
	#
	# Now test whether or not libs for SToGS ROOT Events exists ... to do, call cmake to build it ...
	#
	if (EXISTS ${CMAKE_SOURCE_DIR}/analysis/SToGS/libSToGS_BaseROOTEvents.so)
		#
		ADD_DEFINITIONS(-DHAS_STOGS_ROOT_EVENTS)
		include_directories(${CMAKE_SOURCE_DIR}/analysis/SToGS)
		#
		# corresponding actions managers can be added
    	#set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOTEvent.hh)
		#set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOTEvent.cc)
		#
		#
		set (EXTRA_LIB ${ROOT_LIBRARIES} -L${CMAKE_SOURCE_DIR}/analysis/SToGS -lSToGS_BaseROOTEvents)
		MESSAGE(" --> Cool, libSToGS_BaseROOTEvents.so found in analysis !")	
		#
	else()
		MESSAGE(" --> libSToGS_BaseROOTEvents.so NOT FOUND in analysis, you may needed to compile it first to fully exploit ROOT facilities")	
	endif()
		
    #set (SToGS_HEADER ${SToGS_HEADER} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOTEvent.hh)
    #set (SToGS_SRC ${SToGS_SRC} ${PROJECT_SOURCE_DIR}/csrc/SToGS_BaseROOTEvent.cc)
    #
    MESSAGE("[SToGS] SToGS Compiled with ROOT facilities ")
else()
	MESSAGE("[SToGS] --> ROOT has NOT been found")
endif()
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# ADD My Plugins if any
#
if (MY_DIR)
    # ---------- collect all My files ----------
    set(MY_HEADER)
    file(GLOB MY_HEADER ${MY_DIR}/*.hh ${MY_DIR}/*.h)
    set(MY_SRC)
    file(GLOB MY_SRC ${MY_DIR}/*.cc ${MY_DIR}/*.cpp ${MY_DIR}/*.cxx ${MY_DIR}/*.c)
    #
    include_directories(${PROJECT_SOURCE_DIR}/${MY_DIR})
    #
    #MESSAGE("MY_HEADERS" ${MY_HEADER})
    #MESSAGE("MY_SOURCES" ${MY_SRC})
    
    set (SToGS_HEADER ${SToGS_HEADER} ${MY_HEADER})
    set (SToGS_SRC ${SToGS_SRC} ${MY_SRC})
    
    # ---------- MYDetector ----------
    if (MY_DET)
        if ( EXISTS "${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_DET}DetectorConstruction.hh"
                AND EXISTS "${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_DET}DetectorConstruction.cc" )
            #
            ADD_DEFINITIONS(-DHAS_MYDET)
            #
            MESSAGE("[SToGS] --> ADD ${MY_DET}DetectorConstruction ")
        else()
            MESSAGE("[SToGS] --> Cannot find " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_DET}DetectorConstruction.hh)
            MESSAGE("[SToGS] --> Or " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_DET}DetectorConstruction.cc)
        endif()
    #MESSAGE("[SToGS] --> Compile with ${MY_DET} " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_DET}DetectorConstruction)
    endif()

    # ---------- MYActions ----------
    if (MY_ACT)
        if ( EXISTS "${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_ACT}UserAction.hh"
                AND EXISTS "${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_ACT}UserAction.cc" )
            #
            ADD_DEFINITIONS(-DHAS_MYACT)
            #
            MESSAGE("[SToGS] --> ADD ${MY_ACT}UserAction ")
        else()
            MESSAGE("[SToGS] --> Cannot find " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_ACT}UserAction.hh)
            MESSAGE("[SToGS] --> Or " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_ACT}UserAction.cc)
        endif()
    #MESSAGE("[SToGS] --> Compile with ${MY_ACT} " ${PROJECT_SOURCE_DIR}/${MY_DIR}/${MY_ACT})
    endif()
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
target_link_libraries(SToGS_Source ${Geant4_LIBRARIES} ${SToGS_EXTRA_LIB})
target_link_libraries(SToGS_DetectorBuilder ${Geant4_LIBRARIES} ${SToGS_EXTRA_LIB})


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SToGS_Source DESTINATION bin)
install(TARGETS SToGS_DetectorBuilder DESTINATION bin)



