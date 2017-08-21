CMAKE_MINIMUM_REQUIRED (VERSION 3.8)

PROJECT(gPMC)

SET(GPMC_VERSION_MAJOR "2")
SET(GPMC_VERSION_MINOR "0")
SET(GPMC_VERSION_PATCH "1")

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_AUTOMOC ON)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Wall /favor:INTEL64 /arch:AVX2 /openmp ")

SET(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} /Ox /Wall /favor:INTEL64 /arch:AVX2 /openmp ")

FIND_PACKAGE(OpenCL REQUIRED)
FIND_PACKAGE(Gengetopt)
FIND_PACKAGE(OpenMP)
FIND_PACKAGE(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets DataVisualization REQUIRED)
FIND_PACKAGE(IPP)

IF(OPENMP_FOUND)
  message (STATUS "Compiling with OPENMP option")
  ADD_DEFINITIONS(-DUSE_OPENMP=TRUE)
ENDIF(OPENMP_FOUND)

IF(IPP_FOUND)
  message (STATUS "Compiling with IPP option")
  ADD_DEFINITIONS(-DUSE_IPP=TRUE)
  # SET(IPP_FLAGS)
  include_directories(${IPP_INCLUDE_DIRS})
ENDIF(IPP_FOUND)

SET(ITK_SRC_DIR "${CMAKE_CURRENT_BINARY_DIR}/../gPMCbct/ITKv4" CACHE PATH "ITK source directory")
SET(ITK_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/ITKv120-build" CACHE PATH "ITK build directory")

INCLUDE(ExternalProject)
INCLUDE(cmake/External_ITKv120.cmake)
INCLUDE_DIRECTORIES(${ITK_BUILD_DIR})
INCLUDE_DIRECTORIES(${ITK_BUILD_DIR}/Modules/Core/Common)
FIND_PACKAGE(ITK REQUIRED)
message(STATUS "ITK_DIR = ${ITK_DIR}")
INCLUDE(${ITK_USE_FILE})

SET(GPMC_LIB_DIR "${PROJECT_SOURCE_DIR}/bin")
SET(GPMC_LIBRARIES
  "${GPMC_LIB_DIR}/goPMC.lib"
  "${GPMC_LIB_DIR}/dcmtk.lib"
  "${GPMC_LIB_DIR}/libDicomRT.lib"
)

WRAP_GGO(gPMC_GGO_C gPMC.ggo gPMC3Doutputdose_section.ggo gPMCVersion.ggo.in)

SET(gPMC_FORMS
	gPMC.ui
)
QT5_WRAP_UI(gPMC_FORMS_HEADERS ${gPMC_FORMS})

set(UI_SRC
  gPMC_UI.cxx
  gPMC_UI.h
  AG17RGBAImage.cxx
  AG17RGBAImage.h
  YK16GrayImage.cxx
  YK16GrayImage.h
  qyklabel.cxx
  qyklabel.h
  scatterdatamodifier.cxx
  scatterdatamodifier.h
  range_modulator_data.hxx
  gPMC_commandline.hxx
  gPMC_dcm_tools.hxx
)

ADD_EXECUTABLE(gPMC gPMC.cpp ${gPMC_GGO_C})
ADD_EXECUTABLE(gPMC_UI  main.cpp ${UI_SRC} ${gPMC_FORMS_HEADERS})
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)

MESSAGE(STATUS "CMake of gPMC was found :: ${PROJECT_SOURCE_DIR}")
add_dependencies(gPMC ITKv120)

TARGET_LINK_LIBRARIES(gPMC
  ${GPMC_LIBRARIES}
  ${ITK_LIBRARIES}
  ${IPP_LIBRARIES}
  OpenCL::OpenCL
)

TARGET_LINK_LIBRARIES(gPMC_UI
  ${ITK_LIBRARIES}
  ${IPP_LIBRARIES}
  OpenCL::OpenCL
  Qt5::Core
  Qt5::Gui
  Qt5::OpenGL
  Qt5::PrintSupport
  Qt5::Widgets
  Qt5::DataVisualization
)

SET(GPMC_SRC_DIR ${PROJECT_SOURCE_DIR})

INSTALL( TARGETS gPMC_UI DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin )

SET(GPMC_DLLS
  "${GPMC_LIB_DIR}/goPMC.dll"
  "${GPMC_LIB_DIR}/dcmtk.dll"
  "${GPMC_LIB_DIR}/libDicomRT.dll"
)
INSTALL( FILES ${GPMC_DLLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin )

FILE( GLOB LUT_FILES "${PROJECT_SOURCE_DIR}/input/*.bin" )
INSTALL( FILES ${LUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lut )


if($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets      IMPORTED_LOCATION_DEBUG)
  get_target_property(GUIDLL    Qt5::Gui          IMPORTED_LOCATION_DEBUG)
  get_target_property(COREDLL   Qt5::Core         IMPORTED_LOCATION_DEBUG)
  get_target_property(COREDLL   Qt5::PrintSupport IMPORTED_LOCATION_DEBUG)
  get_target_property(COREDLL   Qt5::DataVisualization IMPORTED_LOCATION_DEBUG)
else($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets      IMPORTED_LOCATION_RELEASE)
  get_target_property(GUIDLL    Qt5::Gui          IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::Core         IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::PrintSupport IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::DataVisualization IMPORTED_LOCATION_RELEASE)
endif($<CONFIG:Debug>)


add_custom_command(
  TARGET gPMC POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:Qt5::Widgets>
    $<TARGET_FILE:Qt5::Gui>
    $<TARGET_FILE:Qt5::Core>
    $<TARGET_FILE:Qt5::PrintSupport>
    $<TARGET_FILE:Qt5::DataVisualization>
    ${CMAKE_CURRENT_BINARY_DIR}/bin
)

add_custom_command(
  TARGET gPMC POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
	${CMAKE_CURRENT_BINARY_DIR}/Debug/gPMC.exe
    ${CMAKE_CURRENT_BINARY_DIR}/bin
)