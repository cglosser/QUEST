find_package (Boost COMPONENTS unit_test_framework)
if(${Boost_FOUND})
  add_executable(qtest-bin "")
  set_target_properties(qtest-bin PROPERTIES OUTPUT_NAME qtest)

  target_compile_features(qtest-bin PUBLIC ${CPP_FEATURES})
  target_sources(qtest-bin PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.cpp")

  include("${CMAKE_CURRENT_LIST_DIR}/AIM/CMakeLists.txt")

  target_sources(qtest-bin
    PUBLIC
      "${CMAKE_CURRENT_LIST_DIR}/integrator_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/lagrange_set_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/math_utils_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/propagator_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/pulse_interaction_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/quantum_dot_test.cpp"
  )

  add_subdirectory("rapidcheck")
  add_subdirectory("rapidcheck/extras/boost_test")

  target_link_libraries(qtest-bin PUBLIC quest ${Boost_LIBRARIES} rapidcheck rapidcheck_boost_test)
endif(${Boost_FOUND})
