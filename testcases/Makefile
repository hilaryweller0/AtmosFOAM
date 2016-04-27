MAKEFLAGS += --no-builtin-rules
.DEFAULT_GOAL := test
.DELETE_ON_ERROR:
.SUFFIXES:
.PHONY: test clean

include $(MAKE_COMMON)/executables/Makefile-OpenFOAM
include $(MAKE_COMMON)/globals/Makefile-OpenFOAM
include $(MAKE_COMMON)/templates/Makefile-FileSystem
include $(MAKE_COMMON)/templates/Makefile-Diagnostics

libfiniteVolumeAtmos := $(FOAM_USER_LIBBIN)/libfiniteVolumeAtmos.so

initialFiles := 0/T 0/phi system/fvSchemes system/fvSolution system/controlDict $(meshFiles)

# $1 -- case directory
# $2 -- l2 error bound
# $3 -- endTime
define AdvectionTest

.PHONY: $1.checkl2errorT

$1.checkl2errorT: build/$1/$3/l2errorBoundT.txt
	@grep -Fxq "1" $$< || (echo "\033[31ml2error exceeded $2\033[0m"; exit 1)

$(eval $(call FieldDiff,build/$1,$3,T,$1/$3/T_analytic))
$(eval $(call L2Error,build/$1,$3,T))
$(eval $(call BoundedL2Error,build/$1,$3,T,$2))

build/$1/$3/T: \
  $(addprefix build/$1/,$(initialFiles)) \
  $(advectionFoam) \
  $(libfiniteVolumeAtmos)
	$(advectionFoam) -case build/$1 > build/$1/advectionFoam.log

$(eval $(call CopyNamedFiles,$1,$(initialFiles),build/$1))

$(eval $(call MakeDir,build/$1/0))
$(eval $(call MakeDir,build/$1/system))
$(eval $(call MakeDir,build/$1/$(polyMesh)))

endef

$(eval $(call AdvectionTest,horizontalBTFCubicUpwindCPCFitAdvection,0.2,5000))
$(eval $(call AdvectionTest,horizontalNoOrographyCubicUpwindCPCFitAdvection,0.006,5000))
$(eval $(call AdvectionTest,terrainFollowingNoOrographyCubicUpwindCPCFitAdvection,0.346,10000))
$(eval $(call AdvectionTest,thermalAdvectionSlantedCellCubicUpwindCPCFit,0.000082,5399.1))
$(eval $(call AdvectionTest,thermalAdvectionSlantedCellFineMeshCubicUpwindCPCFit,0.0000447,3000))
$(eval $(call AdvectionTest,thermalAdvectionBTFCubicUpwindCPCFit,0.00001,5472))
$(eval $(call AdvectionTest,thermalAdvectionCutCellCubicUpwindCPCFit,0.00015,4800))

test: \
  horizontalBTFCubicUpwindCPCFitAdvection.checkl2errorT \
  horizontalNoOrographyCubicUpwindCPCFitAdvection.checkl2errorT \
  terrainFollowingNoOrographyCubicUpwindCPCFitAdvection.checkl2errorT \
  thermalAdvectionSlantedCellCubicUpwindCPCFit.checkl2errorT \
  thermalAdvectionSlantedCellFineMeshCubicUpwindCPCFit.checkl2errorT \
  thermalAdvectionBTFCubicUpwindCPCFit.checkl2errorT \
  thermalAdvectionCutCellCubicUpwindCPCFit.checkl2errorT

clean:
	$(RM) -r build
