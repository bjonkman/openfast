#
# Copyright 2025 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if (GENERATE_TYPES)
  generate_f90_types(src/WaveTank_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/WaveTank_Types.f90 -noextrap)
endif()

#add_library(wavetanktestinglib SHARED 
#    src/WaveTank_Types.f90
#    src/WaveTank_IO.f90
#    src/WaveTank_Struct.f90
#    src/WaveTank.f90
#)
#target_link_libraries(
#    wavetanktestinglib
#        aerodyn_inflow_c_bind_static
#        moordyn_c_bind_static
#        seastate_c_bind_static
#        nwtclibs
#        versioninfolib
#)
#if(APPLE OR UNIX)
#    target_compile_definitions(wavetanktestinglib PRIVATE IMPLICIT_DLLEXPORT)
#endif()
#
#install(TARGETS wavetanktestinglib
#EXPORT "${CMAKE_PROJECT_NAME}Libraries"
#RUNTIME DESTINATION bin
#LIBRARY DESTINATION lib
#ARCHIVE DESTINATION lib)
