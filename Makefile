RM := rm -f
WCLEAN := wclean
WMAKE := FOAM_USER_SRC=${CURDIR} wmake

.PHONY: clean

# TODO: fix dependencies since many executables/libraries depend upon stuff in lnInclude, not just the .so
# we'll probably want .INTERMEDIATE trick again because wmake builds both the .so and lnInclude at the same time

LIB_EXNER_THETA := $(FOAM_USER_LIBBIN)/libExnerTheta.so
LIB_FINITE_VOLUME_USER := $(FOAM_USER_LIBBIN)/libfiniteVolumeUser.so
LIB_FV_MESH_WITH_DUAL := $(FOAM_USER_LIBBIN)/libfvMeshWithDual.so
LIB_HOPS := $(FOAM_USER_LIBBIN)/libHops.so

EXNER_FOAM_H := $(FOAM_USER_APPBIN)/exnerFoamH
ADD_2D_MOUNTAIN := $(FOAM_USER_APPBIN)/add2dMountain
GLOBAL_SUM := $(FOAM_USER_APPBIN)/globalSum
SUM_FIELDS := $(FOAM_USER_APPBIN)/sumFields
SET_EXNER_BALANCED_H := $(FOAM_USER_APPBIN)/setExnerBalancedH
SET_SCALAR_OVER_OROGRAPHY := $(FOAM_USER_APPBIN)/setScalarOverOrography
SET_THETA := $(FOAM_USER_APPBIN)/setTheta
CREATE_SPONGE_LAYER := $(FOAM_USER_APPBIN)/createSpongeLayer
PLOT_PATCH_DATA := $(FOAM_USER_APPBIN)/plotPatchData

ALL_LIBS := $(LIB_EXNER_THETA) \
	$(LIB_FINITE_VOLUME_USER) \
	$(LIB_FV_MESH_WITH_DUAL) \
	$(LIB_HOPS)

ALL_EXECUTABLES := $(EXNER_FOAM_H) \
	$(ADD_2D_MOUNTAIN) \
	$(GLOBAL_SUM) \
	$(SUM_FIELDS) \
	$(SET_EXNER_BALANCED_H) \
	$(SET_SCALAR_OVER_OROGRAPHY) \
	$(SET_THETA) \
	$(CREATE_SPONGE_LAYER) \
	$(PLOT_PATCH_DATA)

all: $(ALL_EXECUTABLES)

clean:
	$(WCLEAN) ExnerTheta
	$(WCLEAN) finiteVolume
	$(WCLEAN) fvMeshWithDual
	$(WCLEAN) Hops
	$(WCLEAN) solvers/ExnerFoam/ExnerFoamH
	$(WCLEAN) solvers/ExnerFoam/ExnerFoam
	$(WCLEAN) utilities/mesh/add2dMountain
	$(WCLEAN) utilities/preProcessing/createSpongeLayer
	$(WCLEAN) utilities/preProcessing/setExnerBalancedH
	$(WCLEAN) utilities/preProcessing/setScalarOverOrography
	$(WCLEAN) utilities/preProcessing/setTheta
	$(WCLEAN) utilities/postProcessing/globalSum
	$(WCLEAN) utilities/postProcessing/sumFields
	$(WCLEAN) utilities/postProcessing/plotPatchData_for5.1.1
	$(RM) $(ALL_LIBS) $(ALL_EXECUTABLES)

$(LIB_EXNER_THETA): ExnerTheta/BCs/fixedFluxBuoyantExnerFvPatchScalarField.C $(LIB_FINITE_VOLUME_USER)
	$(WMAKE) ExnerTheta

# TODO: get dependencies from Make/files ?
$(LIB_FINITE_VOLUME_USER): $(LIB_FV_MESH_WITH_DUAL)
	$(WMAKE) finiteVolume

$(LIB_FV_MESH_WITH_DUAL):
	$(WMAKE) fvMeshWithDual

$(LIB_HOPS):
	$(WMAKE) Hops

$(EXNER_FOAM_H): $(LIB_HOPS) $(LIB_EXNER_THETA)
	$(WMAKE) solvers/ExnerFoam/ExnerFoamH	

$(SET_EXNER_BALANCED_H): $(LIB_HOPS) $(LIB_EXNER_THETA)
	$(WMAKE) utilities/preProcessing/setExnerBalancedH

$(SET_THETA): $(LIB_EXNER_THETA)
	$(WMAKE) utilities/preProcessing/setTheta

$(ADD_2D_MOUNTAIN):
	$(WMAKE) utilities/mesh/add2dMountain

$(GLOBAL_SUM): $(LIB_FV_MESH_WITH_DUAL) $(LIB_FINITE_VOLUME_USER)
	$(WMAKE) utilities/postProcessing/globalSum

$(SUM_FIELDS):
	$(WMAKE) utilities/postProcessing/sumFields

$(SET_SCALAR_OVER_OROGRAPHY):
	$(WMAKE) utilities/preProcessing/setScalarOverOrography

$(CREATE_SPONGE_LAYER): $(LIB_FINITE_VOLUME_USER)
	$(WMAKE) utilities/preProcessing/createSpongeLayer

$(PLOT_PATCH_DATA): $(LIB_FINITE_VOLUME_USER) $(LIB_FV_MESH_WITH_DUAL)
	$(WMAKE) utilities/postProcessing/plotPatchData_for5.1.1
