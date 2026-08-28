cmake_minimum_required(VERSION 3.13)

project(MQMC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(PFQMC_SCALE_SAFE_UDT "Enable scale-safe exponent-aware UDT" ON)
option(PFQMC_BUILD_TESTS "Build the historical GTest suite" OFF)

include_directories(
${PROJECT_SOURCE_DIR}/inc
)
INCLUDE_DIRECTORIES ( ${EIGEN3_INCLUDE_DIR} )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SRC_LIST src/skewMatUtils.cpp src/pfqmc.cpp)
set(PFAPACK_ROOT "${PROJECT_SOURCE_DIR}/inc/pfapack" CACHE PATH
    "Directory containing PFAPACK c_interface/ and fortran/ libraries")

set(MKL_INTERFACE lp64)

set(MKL_ROOT $ENV{MKLROOT})
enable_testing()
find_package(MKL CONFIG QUIET)
if(PFQMC_BUILD_TESTS)
  find_package(GTest CONFIG REQUIRED)
endif()
#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

# add_library( libyaml STATIC IMPORTED )

# set_target_properties( libyaml PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/lib/libyaml-cpp.a )

set(FP_MODEL_FLAGS "-fp-model=precise -fimf-arch-consistency=true -no-fma")

set(WARNING_FLAGS "-O3")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAGS} ${WARNING_FLAGS}")

add_executable(main ./main.cpp ${SRC_LIST})
target_compile_features(main PRIVATE cxx_std_17)
if(PFQMC_SCALE_SAFE_UDT)
  target_compile_definitions(main PRIVATE PFQMC_SCALE_SAFE_UDT)
endif()
if(MKL_FOUND)
  target_compile_options(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_include_directories(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(main PUBLIC $<LINK_ONLY:MKL::MKL>)
else()
  message(STATUS "MKL CMake package unavailable; using Intel -mkl fallback")
  target_compile_options(main PRIVATE -mkl)
  target_include_directories(main PRIVATE $ENV{MKLROOT}/include)
  target_link_libraries(main PUBLIC -mkl)
endif()
target_link_libraries(main PUBLIC ${PFAPACK_ROOT}/c_interface/libcpfapack.a ${PFAPACK_ROOT}/fortran/libpfapack.a)

if(PFQMC_BUILD_TESTS)
  set(TEST_LIST ./test/main_test.cpp ./test/squareLatticeTest.cpp)
  add_executable(main_test ${TEST_LIST} ${SRC_LIST})
  target_compile_features(main_test PRIVATE cxx_std_17)
  if(PFQMC_SCALE_SAFE_UDT)
    target_compile_definitions(main_test PRIVATE PFQMC_SCALE_SAFE_UDT)
  endif()
  if(MKL_FOUND)
    target_compile_options(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(main_test PUBLIC $<LINK_ONLY:MKL::MKL> GTest::gtest_main)
  else()
    target_compile_options(main_test PRIVATE -mkl)
    target_include_directories(main_test PRIVATE $ENV{MKLROOT}/include)
    target_link_libraries(main_test PUBLIC -mkl GTest::gtest_main)
  endif()
  target_link_libraries(main_test PUBLIC ${PFAPACK_ROOT}/c_interface/libcpfapack.a ${PFAPACK_ROOT}/fortran/libpfapack.a)

  include(GoogleTest)
  gtest_discover_tests(main_test)
endif()

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)
