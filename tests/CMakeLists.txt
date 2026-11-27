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
    # Tests tagged [slow] (paper-size inputs) get the label <category>-slow so that
    # `ctest -LE slow` runs the fast suite; CI runs the slow labels in their own step.
    # The split is made only for categories whose sources contain a [slow] tag:
    # CatchAddTests writes no test file when a spec matches nothing, which breaks ctest.
    # A [slow] tag added later without reconfiguring stays under <category> (never dropped).
    file(GLOB_RECURSE _scan CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${category}/*.cpp
         ${CMAKE_CURRENT_SOURCE_DIR}/${category}/*.hpp)
    set(_has_slow FALSE)
    foreach(_f IN LISTS _scan)
      file(STRINGS ${_f} _hit REGEX "\\[slow\\]")
      if(_hit)
        set(_has_slow TRUE)
        break()
      endif()
    endforeach()
    if(_has_slow)
      catch_discover_tests(${_exe} TEST_PREFIX "${category}:" TEST_SPEC "~[slow]" PROPERTIES LABELS ${category} DISCOVERY_MODE PRE_TEST)
      catch_discover_tests(${_exe} TEST_PREFIX "${category}:" TEST_SPEC "[slow]" PROPERTIES LABELS ${category}-slow DISCOVERY_MODE PRE_TEST)
    else()
      catch_discover_tests(${_exe} TEST_PREFIX "${category}:" PROPERTIES LABELS ${category} DISCOVERY_MODE PRE_TEST)
    endif()
  endif()
endforeach()
