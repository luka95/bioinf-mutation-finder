# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\Users\Sime\Develop\CLion\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\Users\Sime\Develop\CLion\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/align.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/align.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/align.dir/flags.make

CMakeFiles/align.dir/Align.cpp.obj: CMakeFiles/align.dir/flags.make
CMakeFiles/align.dir/Align.cpp.obj: ../Align.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/align.dir/Align.cpp.obj"
	D:\Users\Sime\Develop\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\align.dir\Align.cpp.obj -c D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\Align.cpp

CMakeFiles/align.dir/Align.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/align.dir/Align.cpp.i"
	D:\Users\Sime\Develop\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\Align.cpp > CMakeFiles\align.dir\Align.cpp.i

CMakeFiles/align.dir/Align.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/align.dir/Align.cpp.s"
	D:\Users\Sime\Develop\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\Align.cpp -o CMakeFiles\align.dir\Align.cpp.s

# Object files for target align
align_OBJECTS = \
"CMakeFiles/align.dir/Align.cpp.obj"

# External object files for target align
align_EXTERNAL_OBJECTS =

align.exe: CMakeFiles/align.dir/Align.cpp.obj
align.exe: CMakeFiles/align.dir/build.make
align.exe: CMakeFiles/align.dir/linklibs.rsp
align.exe: CMakeFiles/align.dir/objects1.rsp
align.exe: CMakeFiles/align.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable align.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\align.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/align.dir/build: align.exe

.PHONY : CMakeFiles/align.dir/build

CMakeFiles/align.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\align.dir\cmake_clean.cmake
.PHONY : CMakeFiles/align.dir/clean

CMakeFiles/align.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug D:\Users\Sime\FER\Diplomski\Semestar3\Bioinformatika\Projekt\git\bioinf-mutation-finder\cmake-build-debug\CMakeFiles\align.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/align.dir/depend
