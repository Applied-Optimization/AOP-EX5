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
include ConvexityTests/CMakeFiles/ConvexityTests-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include ConvexityTests/CMakeFiles/ConvexityTests-test.dir/compiler_depend.make

# Include the progress variables for this target.
include ConvexityTests/CMakeFiles/ConvexityTests-test.dir/progress.make

# Include the compile flags for this target's objects.
include ConvexityTests/CMakeFiles/ConvexityTests-test.dir/flags.make

ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o: ConvexityTests/CMakeFiles/ConvexityTests-test.dir/flags.make
ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o: /home/deepshukla/aopt-exercise5/aopt-exercise5/ConvexityTests/unit_tests.cc
ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o: ConvexityTests/CMakeFiles/ConvexityTests-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/deepshukla/aopt-exercise5/aopt-exercise5/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o -MF CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o.d -o CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o -c /home/deepshukla/aopt-exercise5/aopt-exercise5/ConvexityTests/unit_tests.cc

ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.i"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/deepshukla/aopt-exercise5/aopt-exercise5/ConvexityTests/unit_tests.cc > CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.i

ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.s"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/deepshukla/aopt-exercise5/aopt-exercise5/ConvexityTests/unit_tests.cc -o CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.s

# Object files for target ConvexityTests-test
ConvexityTests__test_OBJECTS = \
"CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o"

# External object files for target ConvexityTests-test
ConvexityTests__test_EXTERNAL_OBJECTS =

Build/bin/ConvexityTests-test: ConvexityTests/CMakeFiles/ConvexityTests-test.dir/unit_tests.cc.o
Build/bin/ConvexityTests-test: ConvexityTests/CMakeFiles/ConvexityTests-test.dir/build.make
Build/bin/ConvexityTests-test: lib/libgtestd.a
Build/bin/ConvexityTests-test: lib/libgtest_maind.a
Build/bin/ConvexityTests-test: lib/libgtestd.a
Build/bin/ConvexityTests-test: ConvexityTests/CMakeFiles/ConvexityTests-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/deepshukla/aopt-exercise5/aopt-exercise5/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/ConvexityTests-test"
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ConvexityTests-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ConvexityTests/CMakeFiles/ConvexityTests-test.dir/build: Build/bin/ConvexityTests-test
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests-test.dir/build

ConvexityTests/CMakeFiles/ConvexityTests-test.dir/clean:
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests && $(CMAKE_COMMAND) -P CMakeFiles/ConvexityTests-test.dir/cmake_clean.cmake
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests-test.dir/clean

ConvexityTests/CMakeFiles/ConvexityTests-test.dir/depend:
	cd /home/deepshukla/aopt-exercise5/aopt-exercise5/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/deepshukla/aopt-exercise5/aopt-exercise5 /home/deepshukla/aopt-exercise5/aopt-exercise5/ConvexityTests /home/deepshukla/aopt-exercise5/aopt-exercise5/build /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests /home/deepshukla/aopt-exercise5/aopt-exercise5/build/ConvexityTests/CMakeFiles/ConvexityTests-test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests-test.dir/depend

