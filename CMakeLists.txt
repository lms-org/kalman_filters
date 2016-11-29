set(SOURCES
 src/foo.cpp
)

set(HEADERS
    include/kalman_filter/ctra_measurement_model.h
    include/kalman_filter/ctra_system_model.h
    include/kalman_filter/line_ekf.h
    include/kalman_filter/ctra_ekf.h
    include/kalman_filter/mass_2d_rot.h
    include/filter/line_line_x.h
    include/filter/sgd.h
    include/filter/adam.h
    include/utils.h
)

include_directories(
        include
        /usr/local/include
        /usr/include
        /usr/local/include/eigen3
        /usr/include/eigen3
        /usr/include/qt5

)

add_definitions(-DEIGEN_NO_DEBUG)

include_directories(include)
add_library(kalman_filter SHARED ${SOURCES} ${HEADERS})


#TODO just for testing
#QT
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Find the QtWidgets library
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Charts)

add_executable(km
    ${SOURCES} ${HEADERS}
)

target_link_libraries(km PRIVATE lmscore Qt5::Widgets Qt5::Charts Qt5Core)

#target_link_libraries(kalman_filter PRIVATE lmscore)
#add_library(kalman_filter SHARED ${SOURCES} ${HEADERS})
#target_link_libraries(kalman_filter PRIVATE lmscore)
