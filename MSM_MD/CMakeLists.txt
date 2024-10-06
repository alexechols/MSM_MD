# CMakeList.txt : CMake project for MSM_MD, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (MSM_MD "MSM_MD.cpp" "MSM_MD.h" "atoms.h"  "atoms.cpp" "logger.cpp" "logger.h" "utils.cpp" "utils.h" "potential.h" "potential.cpp" "sim.h" "sim.cpp" "input.h" "input.cpp"  "random.h" "random.cpp" )

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET MSM_MD PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
