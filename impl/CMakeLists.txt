# This is the CMake script for compiling a CGAL application.

project (test_01)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(CGAL)
#include(${CGAL_USE_FILE})

if ( CGAL_FOUND )

  include_directories (BEFORE ../../include ./ src/ )

  file(GLOB_RECURSE EMBEDDING_H
    "src/embedding/*.h"
  )
  file(GLOB_RECURSE EMBEDDING_CPP
    "src/embedding/*.cpp"
  )
  file(GLOB_RECURSE OVERLAY_H
    "src/overlay/*.h"
  )
  file(GLOB_RECURSE OVERLAY_CPP
    "src/overlay/*.cpp"
  )
  file(GLOB_RECURSE COMMON_H
    "src/common/*.h"
  )
  file(GLOB_RECURSE COMMON_CPP
    "src/common/*.cpp"
  )
  file(GLOB_RECURSE MASTER_H
    "src/master/*.h"
  )
  file(GLOB_RECURSE MASTER_CPP
    "src/master/*.cpp"
  )
  file(GLOB_RECURSE INTERPOLATION_H
    "src/interpolation/*.h"
  )
  file(GLOB_RECURSE INTERPOLATION_CPP
    "src/interpolation/*.cpp"
  )
  file(GLOB_RECURSE MEASUREMENT_H
    "src/measurement/*.h"
  )
  file(GLOB_RECURSE MEASUREMENT_CPP
    "src/measurement/*.cpp"
  )
  file(GLOB_RECURSE PIPELINE_H
    "src/pipeline/*.h"
  )
  file(GLOB_RECURSE PIPELINE_CPP
    "src/pipeline/*.cpp"
  )
  file(GLOB_RECURSE FEATURE_ALIGNMENT_H
    "src/feature_alignment/*.h"
  )
  file(GLOB_RECURSE FEATURE_ALIGNMENT_CPP
    "src/feature_alignment/*.cpp"
  )
  file(GLOB_RECURSE AUTO_FEATURES_H
    "src/auto_features/*.h"
  )
  file(GLOB_RECURSE AUTO_FEATURES_CPP
    "src/auto_features/*.cpp"
  )
  file(GLOB_RECURSE DIRECT_OPTIMIZATION_H
    "src/direct_optimization/*.h"
  )
  file(GLOB_RECURSE DIRECT_OPTIMIZATION_CPP
    "src/direct_optimization/*.cpp"
  )
  file(GLOB_RECURSE MEASUREMENT_V_H
    "src/measurement_v/*.h"
  )
  file(GLOB_RECURSE MEASUREMENT_V_CPP
    "src/measurement_v/*.cpp"
  )
  
  
  
  # cpp files
  add_executable (morphing
	src/common/constants.h                  
	src/common/typedefs.h
	
												  src/master/main.cpp
  
	src/embedding/main.h 						  src/embedding/main.cpp 
    src/embedding/geometry_utils.h          	  src/embedding/geometry_utils.cpp
    src/embedding/io.h                      	  src/embedding/io.cpp
	src/embedding/debug_geometry.h          	  src/embedding/debug_geometry.cpp
	src/embedding/approx_bounding_sphere.h  	  src/embedding/approx_bounding_sphere.cpp
	src/embedding/random_tetrahedron.h      	  src/embedding/random_tetrahedron.cpp
	src/embedding/relaxation.h              	  src/embedding/relaxation.cpp
	
	src/overlay/main.h 							  src/overlay/main.cpp 
	src/overlay/debug_utils.h					  src/overlay/debug_utils.cpp
    src/overlay/geometry_utils.h            	  src/overlay/geometry_utils.cpp
    src/overlay/cpp_utils.h                 	  src/overlay/cpp_utils.cpp
    src/overlay/impersonator.h              	  src/overlay/impersonator.cpp
    src/overlay/intersections_new.h         	  src/overlay/intersections_new.cpp
    src/overlay/intersections_deprecated.h 	 	  src/overlay/intersections_deprecated.cpp
    src/overlay/mesh_merger.h              	 	  src/overlay/mesh_merger.cpp
    src/overlay/overlay.h              	 	      src/overlay/overlay.cpp
	src/overlay/test_overlay.h          	      src/overlay/test_overlay.cpp
	
	src/interpolation/main.h					  src/interpolation/main.cpp
	src/interpolation/interpolate.h				  src/interpolation/interpolate.cpp
	
	src/measurement/main.h						  src/measurement/main.cpp
	src/measurement/measurement.h				  src/measurement/measurement.cpp
	
	src/pipeline/main.h							  src/pipeline/main.cpp
	
	src/feature_alignment/main.h				  src/feature_alignment/main.cpp
	src/feature_alignment/feature_alignment.h	  src/feature_alignment/feature_alignment.cpp
	src/feature_alignment/grid_minimization.h	  src/feature_alignment/grid_minimization.cpp
	src/feature_alignment/parse_features.h	      src/feature_alignment/parse_features.cpp
	
	src/auto_features/main.h                      src/auto_features/main.cpp
	src/auto_features/simple_features.h           src/auto_features/simple_features.cpp
	src/auto_features/features_io.h               src/auto_features/features_io.cpp
	src/auto_features/low_res_curvature.h         src/auto_features/low_res_curvature.cpp
	
	src/direct_optimization/main.h                src/direct_optimization/main.cpp
	src/direct_optimization/direct_optimization.h src/direct_optimization/direct_optimization.cpp
	
	src/measurement_v/main.h					  src/measurement_v/main.cpp
	src/measurement_v/measurement_v.h			  src/measurement_v/measurement_v.cpp
  )

  source_group("headers\\embedding"           FILES ${EMBEDDING_H})
  source_group("impl\\embedding"              FILES ${EMBEDDING_CPP})
  source_group("headers\\overlay"             FILES ${OVERLAY_H})
  source_group("impl\\overlay"                FILES ${OVERLAY_CPP})
  source_group("headers\\common"              FILES ${COMMON_H})
  source_group("impl\\common"                 FILES ${COMMON_CPP})
  source_group("headers\\master"              FILES ${MASTER_H})
  source_group("impl\\master"                 FILES ${MASTER_CPP})
  source_group("headers\\interpolation"       FILES ${INTERPOLATION_H})
  source_group("impl\\interpolation"          FILES ${INTERPOLATION_CPP})
  source_group("headers\\measurement"         FILES ${MEASUREMENT_H})
  source_group("impl\\measurement"            FILES ${MEASUREMENT_CPP})
  source_group("headers\\pipeline"            FILES ${PIPELINE_H})
  source_group("impl\\pipeline"               FILES ${PIPELINE_CPP})
  source_group("headers\\feature_alignment"   FILES ${FEATURE_ALIGNMENT_H})
  source_group("impl\\feature_alignment"      FILES ${FEATURE_ALIGNMENT_CPP})
  source_group("headers\\auto_features"       FILES ${AUTO_FEATURES_H})
  source_group("impl\\auto_features"          FILES ${AUTO_FEATURES_CPP})
  source_group("headers\\direct_optimization" FILES ${DIRECT_OPTIMIZATION_H})
  source_group("impl\\direct_optimization"    FILES ${DIRECT_OPTIMIZATION_CPP})
  source_group("headers\\measurement_v"       FILES ${MEASUREMENT_V_H})
  source_group("impl\\measurement_v"          FILES ${MEASUREMENT_V_CPP})
  
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS morphing )

  target_link_libraries( morphing PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(morphing)

endif( CGAL_FOUND )
