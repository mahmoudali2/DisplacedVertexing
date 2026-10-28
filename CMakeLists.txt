cmake_minimum_required(VERSION 3.12)
project(DisplacedVertexing LANGUAGES CXX)

# Displaced-vertex finders (arc pair, 3D helix pair, graph-based).
# Needs Gaudi, k4FWCore, EDM4hep and ROOT.

# Same C++ standard as the Key4hep stack
set(CMAKE_CXX_STANDARD 20 CACHE STRING "")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)

# Find required packages
find_package(ROOT REQUIRED COMPONENTS Core MathCore Matrix Physics)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(Gaudi REQUIRED)

# Define sources: framework-free kernels (src/) + Gaudi algorithms (components/)
set(sources
    src/ArcPairVertexing.cpp
    src/HelixPairVertexing.cpp
    src/GraphVertexKernel.cpp
    components/ArcPairVertexFinder.cpp
    components/HelixPairVertexFinder.cpp
    components/GraphDisplacedVertexFinder.cpp
)

# Create the module
gaudi_add_module(DisplacedVertexing
    SOURCES ${sources}
    LINK
        Gaudi::GaudiKernel
        k4FWCore::k4FWCore
        EDM4HEP::edm4hep
        ROOT::Core
        ROOT::MathCore
        ROOT::Matrix
        ROOT::Physics
)

# Include directories
target_include_directories(DisplacedVertexing PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install the target
install(TARGETS DisplacedVertexing
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
)
