cmake_minimum_required (VERSION 3.22)

project(deexGen)

set(SOURCES Angle.cpp      Fus.cpp           LightP.cpp   Random.cpp      TlArray.cpp
AngleDist.cpp  Gdr.cpp           Mass.cpp     Scission.cpp    TlBarDist.cpp
Chart.cpp      History.cpp       Nucleus.cpp  SigBarDist.cpp  Weight.cpp
Evap.cpp       LevelDensity.cpp  Nuclide.cpp  SigCharged.cpp  Yrast.cpp)

set(HEADERS CAngleDist.h  CFus.h           CLightP.h   CRandom.h    CSigBarDist.h  CWeight.h
CAngle.h      CGdr.h           CMass.h     CRun.h       CSigCharged.h  CYrast.h
CChart.h      CHistory.h       CNucleus.h  CRunThick.h  CTlArray.h     SStoreEvap.h
CEvap.h       CLevelDensity.h  CNuclide.h  CScission.h  CTlBarDist.h)

#next we set the ROOT environment
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})#append $ENV{ROOTSYS} to CMAKE_PREFIX_PATH

find_package(ROOT REQUIRED)

if(ROOT_FOUND)
    message("ROOT FOUND")
    include(${ROOT_USE_FILE})
    message(STATUS "ROOT_USE_FILE: ${ROOT_USE_FILE}")
    include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
    message(STATUS "ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}")
    add_definitions(${ROOT_CXX_FLAGS})
    message(STATUS "ROOT_CXX_FLAGS: ${ROOT_CXX_FLAGS}")
endif()

add_executable(deexGen finalGen.cpp ${SOURCES} ${HEADERS})
target_link_libraries(deexGen ${ROOT_LIBRARIES})
message(STATUS "ROOT_LIBRARIES: ${ROOT_LIBRARIES}")

