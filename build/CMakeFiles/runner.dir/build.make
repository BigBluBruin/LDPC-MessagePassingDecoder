# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/GitRepo/LDPC-MessagePassingDecoder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/GitRepo/LDPC-MessagePassingDecoder/build

# Include any dependencies generated for this target.
include CMakeFiles/runner.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/runner.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/runner.dir/flags.make

CMakeFiles/runner.dir/src/Main.cpp.o: CMakeFiles/runner.dir/flags.make
CMakeFiles/runner.dir/src/Main.cpp.o: ../src/Main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/runner.dir/src/Main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runner.dir/src/Main.cpp.o -c /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Main.cpp

CMakeFiles/runner.dir/src/Main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runner.dir/src/Main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Main.cpp > CMakeFiles/runner.dir/src/Main.cpp.i

CMakeFiles/runner.dir/src/Main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runner.dir/src/Main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Main.cpp -o CMakeFiles/runner.dir/src/Main.cpp.s

CMakeFiles/runner.dir/src/Main.cpp.o.requires:

.PHONY : CMakeFiles/runner.dir/src/Main.cpp.o.requires

CMakeFiles/runner.dir/src/Main.cpp.o.provides: CMakeFiles/runner.dir/src/Main.cpp.o.requires
	$(MAKE) -f CMakeFiles/runner.dir/build.make CMakeFiles/runner.dir/src/Main.cpp.o.provides.build
.PHONY : CMakeFiles/runner.dir/src/Main.cpp.o.provides

CMakeFiles/runner.dir/src/Main.cpp.o.provides.build: CMakeFiles/runner.dir/src/Main.cpp.o


CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o: CMakeFiles/runner.dir/flags.make
CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o: ../src/Parity_Check_Matrix_Info.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o -c /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Parity_Check_Matrix_Info.cpp

CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Parity_Check_Matrix_Info.cpp > CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.i

CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/Parity_Check_Matrix_Info.cpp -o CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.s

CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.requires:

.PHONY : CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.requires

CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.provides: CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.requires
	$(MAKE) -f CMakeFiles/runner.dir/build.make CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.provides.build
.PHONY : CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.provides

CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.provides.build: CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o


CMakeFiles/runner.dir/src/full_BP.cpp.o: CMakeFiles/runner.dir/flags.make
CMakeFiles/runner.dir/src/full_BP.cpp.o: ../src/full_BP.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/runner.dir/src/full_BP.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runner.dir/src/full_BP.cpp.o -c /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/full_BP.cpp

CMakeFiles/runner.dir/src/full_BP.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runner.dir/src/full_BP.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/full_BP.cpp > CMakeFiles/runner.dir/src/full_BP.cpp.i

CMakeFiles/runner.dir/src/full_BP.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runner.dir/src/full_BP.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/full_BP.cpp -o CMakeFiles/runner.dir/src/full_BP.cpp.s

CMakeFiles/runner.dir/src/full_BP.cpp.o.requires:

.PHONY : CMakeFiles/runner.dir/src/full_BP.cpp.o.requires

CMakeFiles/runner.dir/src/full_BP.cpp.o.provides: CMakeFiles/runner.dir/src/full_BP.cpp.o.requires
	$(MAKE) -f CMakeFiles/runner.dir/build.make CMakeFiles/runner.dir/src/full_BP.cpp.o.provides.build
.PHONY : CMakeFiles/runner.dir/src/full_BP.cpp.o.provides

CMakeFiles/runner.dir/src/full_BP.cpp.o.provides.build: CMakeFiles/runner.dir/src/full_BP.cpp.o


CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o: CMakeFiles/runner.dir/flags.make
CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o: ../src/ldpc_full_precision_decoder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o -c /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/ldpc_full_precision_decoder.cpp

CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/ldpc_full_precision_decoder.cpp > CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.i

CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/GitRepo/LDPC-MessagePassingDecoder/src/ldpc_full_precision_decoder.cpp -o CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.s

CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.requires:

.PHONY : CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.requires

CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.provides: CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.requires
	$(MAKE) -f CMakeFiles/runner.dir/build.make CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.provides.build
.PHONY : CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.provides

CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.provides.build: CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o


# Object files for target runner
runner_OBJECTS = \
"CMakeFiles/runner.dir/src/Main.cpp.o" \
"CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o" \
"CMakeFiles/runner.dir/src/full_BP.cpp.o" \
"CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o"

# External object files for target runner
runner_EXTERNAL_OBJECTS =

runner: CMakeFiles/runner.dir/src/Main.cpp.o
runner: CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o
runner: CMakeFiles/runner.dir/src/full_BP.cpp.o
runner: CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o
runner: CMakeFiles/runner.dir/build.make
runner: CMakeFiles/runner.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable runner"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/runner.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/runner.dir/build: runner

.PHONY : CMakeFiles/runner.dir/build

CMakeFiles/runner.dir/requires: CMakeFiles/runner.dir/src/Main.cpp.o.requires
CMakeFiles/runner.dir/requires: CMakeFiles/runner.dir/src/Parity_Check_Matrix_Info.cpp.o.requires
CMakeFiles/runner.dir/requires: CMakeFiles/runner.dir/src/full_BP.cpp.o.requires
CMakeFiles/runner.dir/requires: CMakeFiles/runner.dir/src/ldpc_full_precision_decoder.cpp.o.requires

.PHONY : CMakeFiles/runner.dir/requires

CMakeFiles/runner.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/runner.dir/cmake_clean.cmake
.PHONY : CMakeFiles/runner.dir/clean

CMakeFiles/runner.dir/depend:
	cd /mnt/d/GitRepo/LDPC-MessagePassingDecoder/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/GitRepo/LDPC-MessagePassingDecoder /mnt/d/GitRepo/LDPC-MessagePassingDecoder /mnt/d/GitRepo/LDPC-MessagePassingDecoder/build /mnt/d/GitRepo/LDPC-MessagePassingDecoder/build /mnt/d/GitRepo/LDPC-MessagePassingDecoder/build/CMakeFiles/runner.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/runner.dir/depend
