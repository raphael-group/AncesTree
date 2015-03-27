cmake_minimum_required( VERSION 2.8 )

project( AncesTree )

set( analyse_solution_prob_hdr
	src/utils.h
	src/matrix.h
	src/ppmatrix.h
	src/realmatrix.h
	src/maxsolution.h
	src/solutiongraph.h
	src/clonaltree.h
	src/ancestrymatrix.h
	src/baseancestrygraph.h
	src/probancestrygraph.h
	src/readcountmatrix.h
	src/realintervalmatrix.h
	src/probcomparison.h
)

set( analyse_solution_prob_src
	src/analysesolutionprob.cpp
	src/utils.cpp
	src/matrix.cpp
	src/ppmatrix.cpp
	src/realmatrix.cpp
	src/maxsolution.cpp
	src/solutiongraph.cpp
	src/clonaltree.cpp
	src/ancestrymatrix.cpp
	src/baseancestrygraph.cpp
	src/probancestrygraph.cpp
	src/readcountmatrix.cpp
	src/realintervalmatrix.cpp
	src/probcomparison.cpp
)

set( analyse_solution_hdr
	src/utils.h
	src/matrix.h
	src/ppmatrix.h
	src/realmatrix.h
	src/maxsolution.h
	src/solutiongraph.h
	src/clonaltree.h
	src/ancestrymatrix.h
	src/baseancestrygraph.h
	src/probancestrygraph.h
	src/readcountmatrix.h
	src/realintervalmatrix.h
	src/comparison.h
)

set( analyse_solution_src
	src/analysesolution.cpp
	src/utils.cpp
	src/matrix.cpp
	src/ppmatrix.cpp
	src/realmatrix.cpp
	src/maxsolution.cpp
	src/solutiongraph.cpp
	src/clonaltree.cpp
	src/ancestrymatrix.cpp
	src/baseancestrygraph.cpp
	src/probancestrygraph.cpp
	src/readcountmatrix.cpp
	src/realintervalmatrix.cpp
	src/comparison.cpp
)

set( visualize_solution_src
	src/visualizesolution.cpp
	src/utils.cpp
	src/matrix.cpp
	src/ppmatrix.cpp
	src/realmatrix.cpp
	src/maxsolution.cpp
	src/solutiongraph.cpp
	src/clonaltree.cpp
)

set( visualize_solution_hdr
	src/utils.h
	src/matrix.h
	src/ppmatrix.h
	src/realmatrix.h
	src/maxsolution.h
	src/solutiongraph.h
	src/clonaltree.h
)

set( analyse_src
	src/analyse.cpp
	src/probancestrygraph.cpp
	src/baseancestrygraph.cpp
	src/ancestrymatrix.cpp
	src/utils.cpp
	src/matrix.cpp
	src/readcountmatrix.cpp
	src/realintervalmatrix.cpp
	src/realmatrix.cpp
)

set( analyse_hdr
	src/probancestrygraph.h
	src/baseancestrygraph.h
	src/ancestrymatrix.h
	src/utils.h
	src/matrix.h
	src/readcountmatrix.h
	src/realintervalmatrix.h
	src/realmatrix.h
)

set( ancestree_ilp_src
	src/ancestreeilp.cpp
	src/intmaxilpsolver.cpp
	src/probancestrygraph.cpp
	src/baseancestrygraph.cpp
	src/ancestrymatrix.cpp
	src/utils.cpp
	src/matrix.cpp
	src/ppmatrix.cpp
	src/realmatrix.cpp
	src/realintervalmatrix.cpp
	src/readcountmatrix.cpp
	src/maxsolution.cpp
	src/clonaltree.cpp
)

set( ancestree_ilp_hdr
	src/intmaxilpsolver.h
	src/probancestrygraph.h
	src/ancestrymatrix.h
	src/baseancestrygraph.h
	src/utils.h
	src/matrix.h
	src/ppmatrix.h
	src/realmatrix.h
	src/realintervalmatrix.h
	src/readcountmatrix.h
	src/maxsolution.h
	src/clonaltree.h
)

set( constructancestrymatrix_src
	src/constructancestrymatrix.cpp
	src/utils.cpp
	src/readcountmatrix.cpp
	src/matrix.cpp
	src/ancestrymatrix.cpp
	src/realintervalmatrix.cpp
	src/realmatrix.cpp
)

set( constructancestrymatrix_hdr
	src/utils.h
	src/readcountmatrix.h
	src/matrix.h
	src/realintervalmatrix.h
	src/realmatrix.h
	src/ancestrymatrix.h
)

set( ancestree_src
	src/ancestree.cpp
	src/solutiongraph.cpp
	src/intmaxilpsolver.cpp
	src/probancestrygraph.cpp
	src/baseancestrygraph.cpp
	src/ppmatrix.cpp
	src/maxsolution.cpp
	src/clonaltree.cpp
	src/utils.cpp
	src/readcountmatrix.cpp
	src/matrix.cpp
	src/realintervalmatrix.cpp
	src/realmatrix.cpp
	src/ancestrymatrix.cpp
)

set( ancestree_hdr
	src/solutiongraph.h
	src/intmaxilpsolver.h
	src/probancestrygraph.h
	src/baseancestrygraph.h
	src/ppmatrix.h
	src/maxsolution.h
	src/clonaltree.h
	src/utils.h
	src/readcountmatrix.h
	src/matrix.h
	src/realintervalmatrix.h
	src/realmatrix.h
	src/ancestrymatrix.h
	src/config.h
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( ANCESTREE_VERSION \"${GitBranch}-${GitHashtag}\" )

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost REQUIRED )

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
	set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
	foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
		file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
		if( CplexLibraryPaths )
			list( LENGTH CplexLibraryPaths NumPath )
			if( NumPath GREATER 1 )
				message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
			else()
				message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
			endif()
			list( GET CplexLibraryPaths 0 CplexLibraryPath )

			set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
			string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
			if( NOT ValidatedCplexLibPath )
				message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
			else()
				string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
				string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

				set( CPLEX_INC_DIR	 "${CplexInstallPrefix}/${CplexHome}/cplex/include"						CACHE PATH "CPLEX's include dir" )
				set( CPLEX_LIB_DIR	 "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"	 CACHE PATH "CPLEX's library dir" )
				set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"					CACHE PATH "Concert's include dir" )
				set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

				break()
			endif()
		else()
			message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
		endif()
	endforeach()

	if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
		message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
	else()
		message( "Autodetected CPLEX's include dir	 : ${CPLEX_INC_DIR}" )
		message( "Autodetected CPLEX's library dir	 : ${CPLEX_LIB_DIR}" )
		message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
		message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
	endif()
endif()

SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
add_definitions( -DIL_STD )
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
	m
	ilocplex
	cplex
	concert
	pthread
	emon
)

add_executable( analyse_solution ${analyse_solution_src} ${analyse_solution_hdr} )
target_link_libraries( analyse_solution ${CommonLibs} )

add_executable( analyse_solution_prob ${analyse_solution_prob_src} ${analyse_solution_prob_hdr} )
target_link_libraries( analyse_solution_prob ${CommonLibs} )

add_executable( visualize_solution ${visualize_solution_src} ${visualize_solution_hdr} )
target_link_libraries( visualize_solution ${CommonLibs} )

add_executable( construct_ancestry_matrix ${constructancestrymatrix_src} ${constructancestrymatrix_hdr} )
target_link_libraries( construct_ancestry_matrix ${CommonLibs} )

add_executable( analyse EXCLUDE_FROM_ALL ${analyse_src} ${analyse_hdr} )
target_link_libraries( analyse ${CommonLibs} )

add_executable( ancestree_ilp ${ancestree_ilp_src} ${ancestree_ilp_hdr} )
target_link_libraries( ancestree_ilp ${CommonLibs} )

add_executable( ancestree ${ancestree_src} ${ancestree_hdr} )
target_link_libraries( ancestree ${CommonLibs} )
