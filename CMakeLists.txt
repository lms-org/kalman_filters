set(SOURCES
 src/foo.cpp
)

set(HEADERS
    include/kalman_filter/ctra_measurement_model.h
    include/kalman_filter/ctra_system_model.h
    include/kalman_filter/line_ekf.h
    include/kalman_filter/ctra_ekf.h
    include/kalman_filter/mass_ekf.h
    include/filter/line_line_x.h
    include/utils.h
)

include_directories(
        include
        /usr/local/include
        /usr/include
        /usr/local/include/eigen3
        /usr/include/eigen3
)

add_definitions(-DEIGEN_NO_DEBUG)

include_directories(include)
add_library(kalman_filter SHARED ${SOURCES} ${HEADERS})
#target_link_libraries(kalman_filter PRIVATE lmscore)
