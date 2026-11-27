include(Catch)
# One executable per test category; files are globbed so modules can be added
# in parallel without touching this file.
foreach(category unit integration differential property regression)
  file(GLOB_RECURSE _srcs CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${category}/*.cpp)
  if(_srcs)
    set(_exe sssp_${category}_tests)
    add_executable(${_exe} ${_srcs})
    target_link_libraries(${_exe} PRIVATE sssp Catch2::Catch2WithMain)
    target_include_directories(${_exe} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src)
    sssp_set_warnings(${_exe})
    catch_discover_tests(${_exe} TEST_PREFIX "${category}:" PROPERTIES LABELS ${category} DISCOVERY_MODE PRE_TEST)
  endif()
endforeach()
