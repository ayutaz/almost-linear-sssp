function(sssp_set_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
  if(SSSP_WERROR)
    target_compile_options(${target} PRIVATE -Werror)
  endif()
  # fast-math is never enabled: correctness relies on IEEE-754 semantics.
  target_compile_options(${target} PRIVATE -fno-fast-math)
endfunction()

if(SSSP_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined)
  add_link_options(-fsanitize=address,undefined)
endif()

if(SSSP_COVERAGE)
  add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
  add_link_options(-fprofile-instr-generate -fcoverage-mapping)
endif()
