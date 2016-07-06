set(SOURCES
# src/file.cpp
)

set(HEADERS
# include/kalman_filter/file.h
)

include_directories(include)
add_library(kalman_filter SHARED ${SOURCES} ${HEADERS})
target_link_libraries(kalman_filter PRIVATE lmscore)
