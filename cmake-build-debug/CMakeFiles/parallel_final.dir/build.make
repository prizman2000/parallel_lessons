# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/apple/CLionProjects/parallel_final

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/apple/CLionProjects/parallel_final/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/parallel_final.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/parallel_final.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/parallel_final.dir/flags.make

CMakeFiles/parallel_final.dir/main.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/parallel_final.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/main.cpp.o -c /Users/apple/CLionProjects/parallel_final/main.cpp

CMakeFiles/parallel_final.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/main.cpp > CMakeFiles/parallel_final.dir/main.cpp.i

CMakeFiles/parallel_final.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/main.cpp -o CMakeFiles/parallel_final.dir/main.cpp.s

CMakeFiles/parallel_final.dir/src/threads_num.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/src/threads_num.cpp.o: ../src/threads_num.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/parallel_final.dir/src/threads_num.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/src/threads_num.cpp.o -c /Users/apple/CLionProjects/parallel_final/src/threads_num.cpp

CMakeFiles/parallel_final.dir/src/threads_num.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/src/threads_num.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/src/threads_num.cpp > CMakeFiles/parallel_final.dir/src/threads_num.cpp.i

CMakeFiles/parallel_final.dir/src/threads_num.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/src/threads_num.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/src/threads_num.cpp -o CMakeFiles/parallel_final.dir/src/threads_num.cpp.s

CMakeFiles/parallel_final.dir/src/functions.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/src/functions.cpp.o: ../src/functions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/parallel_final.dir/src/functions.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/src/functions.cpp.o -c /Users/apple/CLionProjects/parallel_final/src/functions.cpp

CMakeFiles/parallel_final.dir/src/functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/src/functions.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/src/functions.cpp > CMakeFiles/parallel_final.dir/src/functions.cpp.i

CMakeFiles/parallel_final.dir/src/functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/src/functions.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/src/functions.cpp -o CMakeFiles/parallel_final.dir/src/functions.cpp.s

CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o: ../src/integrate_cpp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o -c /Users/apple/CLionProjects/parallel_final/src/integrate_cpp.cpp

CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/src/integrate_cpp.cpp > CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.i

CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/src/integrate_cpp.cpp -o CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.s

CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o: ../src/integrate_omp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o -c /Users/apple/CLionProjects/parallel_final/src/integrate_omp.cpp

CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/src/integrate_omp.cpp > CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.i

CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/src/integrate_omp.cpp -o CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.s

CMakeFiles/parallel_final.dir/src/reduction.cpp.o: CMakeFiles/parallel_final.dir/flags.make
CMakeFiles/parallel_final.dir/src/reduction.cpp.o: ../src/reduction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/parallel_final.dir/src/reduction.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel_final.dir/src/reduction.cpp.o -c /Users/apple/CLionProjects/parallel_final/src/reduction.cpp

CMakeFiles/parallel_final.dir/src/reduction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel_final.dir/src/reduction.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/apple/CLionProjects/parallel_final/src/reduction.cpp > CMakeFiles/parallel_final.dir/src/reduction.cpp.i

CMakeFiles/parallel_final.dir/src/reduction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel_final.dir/src/reduction.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/apple/CLionProjects/parallel_final/src/reduction.cpp -o CMakeFiles/parallel_final.dir/src/reduction.cpp.s

# Object files for target parallel_final
parallel_final_OBJECTS = \
"CMakeFiles/parallel_final.dir/main.cpp.o" \
"CMakeFiles/parallel_final.dir/src/threads_num.cpp.o" \
"CMakeFiles/parallel_final.dir/src/functions.cpp.o" \
"CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o" \
"CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o" \
"CMakeFiles/parallel_final.dir/src/reduction.cpp.o"

# External object files for target parallel_final
parallel_final_EXTERNAL_OBJECTS =

parallel_final: CMakeFiles/parallel_final.dir/main.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/src/threads_num.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/src/functions.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/src/integrate_cpp.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/src/integrate_omp.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/src/reduction.cpp.o
parallel_final: CMakeFiles/parallel_final.dir/build.make
parallel_final: /usr/local/lib/libomp.dylib
parallel_final: CMakeFiles/parallel_final.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable parallel_final"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parallel_final.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/parallel_final.dir/build: parallel_final

.PHONY : CMakeFiles/parallel_final.dir/build

CMakeFiles/parallel_final.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/parallel_final.dir/cmake_clean.cmake
.PHONY : CMakeFiles/parallel_final.dir/clean

CMakeFiles/parallel_final.dir/depend:
	cd /Users/apple/CLionProjects/parallel_final/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/CLionProjects/parallel_final /Users/apple/CLionProjects/parallel_final /Users/apple/CLionProjects/parallel_final/cmake-build-debug /Users/apple/CLionProjects/parallel_final/cmake-build-debug /Users/apple/CLionProjects/parallel_final/cmake-build-debug/CMakeFiles/parallel_final.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/parallel_final.dir/depend
