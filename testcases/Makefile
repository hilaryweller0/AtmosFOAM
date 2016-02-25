MAKEFLAGS += --no-builtin-rules
.DEFAULT_GOAL := test
.DELETE_ON_ERROR:
.SUFFIXES:
.PHONY: test

include $(MAKE_COMMON)/globals/Makefile-OpenFOAM
include $(MAKE_COMMON)/templates/Makefile-FileSystem

advectionFoam := $(FOAM_USER_APPBIN)/advectionFoam
libfiniteVolumeAtmos := $(FOAM_USER_LIBBIN)/libfiniteVolumeAtmos.so

initialFiles := 0/T 0/phi system/fvSchemes system/fvSolution system/controlDict $(meshFiles)

test: build/horizontalCubicUpwindCPCFitAdvection/5000/T

build/horizontalCubicUpwindCPCFitAdvection/5000/T: \
  $(addprefix build/horizontalCubicUpwindCPCFitAdvection/,$(initialFiles)) \
  $(advectionFoam) \
  $(libfiniteVolumeAtmos)
	$(advectionFoam) -case build/horizontalCubicUpwindCPCFitAdvection

$(eval $(call CopyNamedFiles,horizontalCubicUpwindCPCFitAdvection,$(initialFiles),build/horizontalCubicUpwindCPCFitAdvection))

$(eval $(call MakeDir,build/horizontalCubicUpwindCPCFitAdvection/0))
$(eval $(call MakeDir,build/horizontalCubicUpwindCPCFitAdvection/system))
$(eval $(call MakeDir,build/horizontalCubicUpwindCPCFitAdvection/$(polyMesh)))
