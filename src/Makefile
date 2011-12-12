# I (Michael Eckmann) did not create this original Makefile
# I got it from somewhere but from where, I don't recall.
# If I do ever run across where I found it, I will make sure I am allowed
# to use it AND give proper credit to the author.
#
# That said, all that needs to be done is to separate the program names
# by spaces in the ALL_PROGRAMS line below and add a line for
# any added programs with the name of the program followed immediately
# by _OBJS = 
# then all the .o files separated by whitespace that are to be created
# and linked together to form the program

# The programs we would like to build
ALL_PROGRAMS = AddImagesTest

# The object files we need to build each of the programs
AddImagesTest_OBJS = AddImagesTest.o

# Location of OPENCV libraries
OPENCV	= /usr/lib

# Location of OPENCV include files
OPENCV_INC	= /usr/include/opencv

# C++ compiler is traditionally called CXX
CXX = g++

# The traditional name for the flags to pass to the preprocessor.
# Useful if you need to specify additional directories to look for
# header files.
#CPPFLAGS = # -I../my_headers
CPPFLAGS	= -g -I$(OPENCV_INC) -Wall -ansi

# A traditional name for the C++ compiler flags
CXXFLAGS = -g -Wall

# Command to make an object file:
COMPILE = $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c

# Command used to link objects:
# LD is sort of traditional as well,
# probably stands for <b>L</b>inker-loa<b>D</b>er
LD = $(CXX)

# The flags to pass to the linker
# At UW CSE, if you are using /uns/g++ as your compiler,
# You'll probably want to add -Wl,-rpath,/uns/lib here.
# This way you won't have to set the LD_LIBRARY_PATH environment
# variable just to run the programs you compile.
#LDFLAGS = # -Wl,-rpath,/uns/lib
LDFLAGS	= -g -L$(OPENCV)/lib

# Any libraries we may want to use
# you would need a -lm here if you were using math functions.
LIBS	= -lm -lcv -lhighgui -lcvaux

# The rule that occurs first in the makefile is the default.
# By default, we want to build everything.
all : $(ALL_PROGRAMS)

# We don't want Make to get confused if a file named 
# "all" should happen to exist, so we say that
# "all" is a PHONY target, i.e., a target that Make
# should always try to update.  We do that by making all a
# dependency of the special target .PHONY:
.PHONY : all


## What is it the user wishes to build?
# The special variable MAKECMDGOALS contains the list of targets the
# user has supplied on the "make" command line.  The
# MAKECMDGOALS variable will only be "defined" (will have a
# value) if the user had actually provided a target on the command
# line.
ifdef MAKECMDGOALS
# the user specified a target on the command line:
# was it "all", or was it something else?
ifeq "$(MAKECMDGOALS)" "all"
programs_to_build = $(ALL_PROGRAMS)
else
programs_to_build = $(MAKECMDGOALS)
endif # ? "$(MAKECMDGOALS)" == "all"
else # ! MAKECMDGOALS
# the user didn't specify a target on the command line 
# (i.e., they just typed "make"):
# act as if they had specified the "all" target
programs_to_build = $(ALL_PROGRAMS)
endif # ? MAKECMDGOALS

# Now that we know what we are building, figure out which .link and .d
# files we need to pull in.
dot_link_files = $(programs_to_build:%=%.link)
objects_to_compile = $(sort $(foreach program,$(programs_to_build),$($(program)_OBJS)))
dot_d_files = $(objects_to_compile:.o=.d)

# Now pull in those .link and .d files -- but only if the list of
# objects we need to compile is not empty.
ifneq "$(strip $(objects_to_compile))" ""
-include $(dot_link_files)
-include $(dot_d_files)
endif

# Here is how to make a .link mini-makefile which gives a rule to link
# a particular program
%.link : always_remake
	echo '$* : $($*_OBJS)' >$@
	echo '#$$(LD) $$(LDFLAGS) $$^ $$(LIBS) -o $$@' \
		| tr '#' '\t' >>$@

# The only file a .link file really properly depends on is the
# makefile itself.  However, I could not find a way to tell
# which makefile GNU make is currently processing.  So, instead,
# I always force the .link file to be remade.  This is a
# drawback of my apporach.
.PHONY : always_remake


# Here is how to produce a .d mini-makefile, listing the 
# dependencies of a particular source file
%.d : %.cpp
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
		| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@

# Here is how we want to compile our source files
%.o : %.cpp
	$(COMPILE) $< -o $@

# A rule to clean up the build directory
.PHONY : clean
clean:
	rm -f *.o *.link *.d $(ALL_PROGRAMS)
