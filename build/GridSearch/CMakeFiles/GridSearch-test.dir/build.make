# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/deepshukla/aopt-exercise5/aopt-exercise5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/deepshukla/aopt-exercise5/aopt-exercise5/build

# Include any dependencies generated for this target.
include GridSearch/CMakeFiles/GridSearch-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include GridSearch/CMakeFiles/GridSearch-test.dir/compiler_depend.make

# Include the progress variables for this target.
include GridSearch/CMakeFiles/GridSearch-test.dir/progress.make

# Include the compile flags for this target's objects.
include GridSearch/CMakeFiles/GridSearch-test.dir/flags.make

GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o: GridSearch/CMakeFiles/GridSearch-test.dir/flags.make
GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o: /home/deepshukla/aopt-exercise5/aopt-exercise5/GridSearch/unit_tests.cc
GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o: GridSearch/CMakeFiles/GridSearch-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/deepshukla/aopt-exercise5/aopt-exercise5/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o -MF CMakeFiles/GridSearch-test.dir/unit_tests.cc.o.d -o CMakeFiles/GridSearch-test.dir/unit_tests.cc.o -c /home/deepshukla/aopt-exercise5/aopt-exercise5/GridSearch/unit_tests.cc

GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/GridSearch-test.dir/unit_tests.cc.i"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/deepshukla/aopt-exercise5/aopt-exercise5/GridSearch/unit_tests.cc > CMakeFiles/GridSearch-test.dir/unit_tests.cc.i

GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/GridSearch-test.dir/unit_tests.cc.s"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/deepshukla/aopt-exercise5/aopt-exercise5/GridSearch/unit_tests.cc -o CMakeFiles/GridSearch-test.dir/unit_tests.cc.s

# Object files for target GridSearch-test
GridSearch__test_OBJECTS = \
"CMakeFiles/GridSearch-test.dir/unit_tests.cc.o"

# External object files for target GridSearch-test
GridSearch__test_EXTERNAL_OBJECTS =

Build/bin/GridSearch-test: GridSearch/CMakeFiles/GridSearch-test.dir/unit_tests.cc.o
Build/bin/GridSearch-test: GridSearch/CMakeFiles/GridSearch-test.dir/build.make
Build/bin/GridSearch-test: lib/libgtestd.a
Build/bin/GridSearch-test: lib/libgtest_maind.a
Build/bin/GridSearch-test: lib/libgtestd.a
Build/bin/GridSearch-test: GridSearch/CMakeFiles/GridSearch-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/deepshukla/aopt-exercise5/aopt-exercise5/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/GridSearch-test"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GridSearch-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
GridSearch/CMakeFiles/GridSearch-test.dir/build: Build/bin/GridSearch-test
.PHONY : GridSearch/CMakeFiles/GridSearch-test.dir/build

GridSearch/CMakeFiles/GridSearch-test.dir/clean:
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch && $(CMAKE_COMMAND) -P CMakeFiles/GridSearch-test.dir/cmake_clean.cmake
.PHONY : GridSearch/CMakeFiles/GridSearch-test.dir/clean

GridSearch/CMakeFiles/GridSearch-test.dir/depend:
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/deepshukla/aopt-exercise5/aopt-exercise5 /home/deepshukla/aopt-exercise5/aopt-exercise5/GridSearch /home/deepshukla/aopt-exercise5/aopt-exercise5/build /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch /home/deepshukla/aopt-exercise5/aopt-exercise5/build/GridSearch/CMakeFiles/GridSearch-test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : GridSearch/CMakeFiles/GridSearch-test.dir/depend

