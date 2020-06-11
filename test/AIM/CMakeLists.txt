target_sources(qtest-bin
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/aim_test.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/expansion_test.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/grid_test.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/spacetime_test.cpp"

    "${CMAKE_CURRENT_LIST_DIR}/farfield_test.cpp"
    #"${CMAKE_CURRENT_LIST_DIR}/field_tests.cpp"
)
