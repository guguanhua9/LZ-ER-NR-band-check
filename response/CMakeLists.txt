cmake_minimum_required(VERSION 3.16)
project(LZBandResponseSnapshot LANGUAGES CXX)

# A local frozen reproduction subset, not NEST's official build definition.
# All includes and translation units are bundled. No FetchContent/network,
# prebuilt .so/.a, ROOT, Geant4, pybind11, or installed NEST is used.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(NEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/nest")
set(GCEM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/gcem")
add_library(NESTCore STATIC
  "${NEST_ROOT}/src/NEST.cpp"
  "${NEST_ROOT}/src/RandomGen.cpp"
  "${NEST_ROOT}/src/TestSpectra.cpp"
  "${NEST_ROOT}/src/VDetector.cpp"
  "${NEST_ROOT}/src/GammaHandler.cpp"
  "${NEST_ROOT}/src/ValidityTests.cpp"
  "${NEST_ROOT}/src/LArNEST.cpp"
  "${NEST_ROOT}/src/LArDetector.cpp"
  "${NEST_ROOT}/src/TestLArSpectra.cpp"
)
target_include_directories(NESTCore PUBLIC
  "${NEST_ROOT}/include"
  "${NEST_ROOT}/include/NEST"
  "${NEST_ROOT}/include/Detectors"
  "${GCEM_ROOT}/include"
)
add_executable(lz_response src/run_response.cpp src/nestDiamxCore.cpp)
target_link_libraries(lz_response PRIVATE NESTCore)
target_include_directories(lz_response PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
