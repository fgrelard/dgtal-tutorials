FILE(REMOVE_RECURSE
  "CMakeFiles/Hello.dir/Hello.cpp.o"
  "Hello.pdb"
  "Hello"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang CXX)
  INCLUDE(CMakeFiles/Hello.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
