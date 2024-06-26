#-----------------------------------------------------------------------
# This Makefile is for building MOZART2 Pre-processor
#------------------------------------------------------------------------

# Set up special characters
null  :=
space := $(null) $(null)

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := ../bin
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := campp
endif

# Check for source list
ifeq ($(SRCLIST),$(null))
SRCLIST := Base_Srclist_f
endif

ifeq ($(OBJ_DIR),$(null))
OBJ_DIR := ./OBJ
endif

# Load dependency search path.
# Check for source directories
ifeq ($(SRCDIRS),$(null))
dirs := ./cam_chempp
else
dirs := ./ $(SRCDIRS)
endif

# Determine platform
UNAMES := $(shell uname -s)
UNAMEM := $(shell uname -m)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(dir))
VPATH    := $(subst $(space),:,$(VPATH))

# Get list of files and determine objects and dependency files
base_srclist_f := $(shell cat $(SRCLIST))
OBJS      := $(foreach file,$(base_srclist_f),$(OBJ_DIR)/$(file:.f=.o))

all: $(MODEL_EXEDIR)/$(EXENAME)

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------

# guess default compiler
ifeq ($(USER_FC),$(null))

#------------------------------------------------------------------------
# Linux -- including pleiades
#------------------------------------------------------------------------
ifeq ($(UNAMES),Linux)
  USER_FC := gfortran
endif
#------------------------------------------------------------------------
# Altix (columbia)
#------------------------------------------------------------------------
ifeq ($(UNAMEM),ia64)
  USER_FC := ifort
endif
#------------------------------------------------------------------------
# Alpha
#------------------------------------------------------------------------
ifeq ($(UNAMEM),alpha)
  USER_FC := f95
endif
#------------------------------------------------------------------------
# AIX ?
#------------------------------------------------------------------------
ifeq ($(UNAMES),AIX)
  USER_FC := xlf95
endif
#------------------------------------------------------------------------
# BG/L, BG/P
#------------------------------------------------------------------------
ifeq ($(UNAMEM),ppc64)
  USER_FC := xlf95
endif
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)
  USER_FC := f90
  MACHFLGS := -OPT:Olimit=8200
endif

#------------------------------------------------------------------------
#------------------------------------------------------------------------

endif

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# set compiler flags ...
#------------------------------------------------------------------------
ifeq ($(USER_FC),ifort)
  FFLAGS      :=  -O2 -c -132 -ftz -g -FR  -I $(OBJ_DIR)
  ifeq ($(DEBUG),TRUE)
    FFLAGS    += -CB
  endif
endif
ifeq ($(USER_FC),lf95)
  ifeq ($(DEBUG),TRUE)
    FFLAGS    := --nfix --nap --chk --g --npca --nsav --trace --trap -c --mod $(OBJ_DIR) -O
  else
    FFLAGS    := --nfix --nap --nchk --ng --npca --nsav --ntrace -c --mod $(OBJ_DIR) -O2
  endif
endif
ifeq ($(USER_FC),pgf90)
  FFLAGS      := -O1 -c -g -Mfree -module $(OBJ_DIR)
  ifeq ($(DEBUG),TRUE)
    FFLAGS    += -C
  endif
endif
ifeq ($(USER_FC),pgf95)
  FFLAGS      := -O1 -c -g -Mfree -module $(OBJ_DIR)
  ifeq ($(DEBUG),TRUE)
    FFLAGS    += -C
  endif
endif
ifeq ($(USER_FC),f90)
  FFLAGS      := -c -freeform -I $(OBJ_DIR) -O2 $(MACHFLGS)
endif
ifeq ($(USER_FC),f95)
  FFLAGS      := -O4 -c -tune host -arch host -free -module $(OBJ_DIR) -I $(OBJ_DIR)
endif
ifeq ($(USER_FC),xlf95)
  FFLAGS      := -g -c -qarch=auto -qnosave -qfree=f90 -qmoddir=$(OBJ_DIR) -I $(OBJ_DIR) -qstrict -O3
endif
ifeq ($(USER_FC),gfortran)
  FFLAGS      := -g -c -ffree-form
endif
ifeq ($(USER_FC),g95)
  FFLAGS      := -g -c -ffree-form
endif

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------

FC := $(USER_FC)

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) -o $@ $<

$(MODEL_EXEDIR)/$(EXENAME):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

RM := rm

clean:
	$(RM) -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(MODEL_EXEDIR)/$(EXENAME)

realclean:
	$(RM) -f $(OBJ_DIR)/*.o *.d $(MODEL_EXEDIR)/$(EXENAME)
