add_executable (main main.cc)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)   # See below (1)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(main GSL::gsl GSL::gslcblas)