# LAYNII makefile

CC		= g++
CFLAGS	= -std=c++11 -DHAVE_ZLIB -lm -lz
# CFLAGS	= -std=c++11 -pedantic -DHAVE_ZLIB -lm -lz

# =============================================================================
LIBRARIES		=	dep/nifti1_io.cpp \
					dep/znzlib.cpp \
					dep/laynii_lib.cpp \
					-I./niftilib \
					-I./znzlib \
					-I./laynii_lib

HIGH_PRIORITY	= 	LN_BOCO \
					LN_MP2RAGE_DNOISE \
					LN2_LAYER_SMOOTH \
					LN_LAYER_SMOOTH \

LOW_PRIORITY	=	LN_3DCOLUMNS \
					LN_COLUMNAR_DIST \
					LN_CORREL2FILES \
					LN_DIRECT_SMOOTH \
					LN_GRADSMOOTH \
					LN_ZOOM \
					LN_FLOAT_ME \
					LN_SHORT_ME \
					LN_EXTREMETR \
					LN_GFACTOR \
					LN_GROW_LAYERS \
					LN_IMAGIRO \
					LN_INTPRO \
					LN_LEAKY_LAYERS \
					LN_NOISEME \
					LN_RAGRUG \
					LN_SKEW \
					LN_TEMPSMOOTH \
					LN_TRIAL \
					LN_PHYSIO_PARS \
					LN_INT_ME \
					LN_LOITUMA \
					LN_NOISE_KERNEL \
					LN2_DEVEIN \

LAYNII2	= LN2_LAYERS

LAYNII 	= $(HIGH_PRIORITY) $(LOW_PRIORITY) $(LAYNII2)

# =============================================================================
all : $(LAYNII)

.PHONY: all $(HIGH_PRIORITY) $(LOW_PRIORITY) $(LAYNII2)

# =============================================================================
# LAYNII v2.0.0 programs
LN2_LAYERS:
	$(CC) $(CFLAGS) -o LN2_LAYERS src/LN2_LAYERS.cpp $(LIBRARIES)

LN2_LLOYD:
	$(CC) $(CFLAGS) -o obj/LN2_LLOYD.o src/LN2_LLOYD.cpp $(LIBRARIES)

# -----------------------------------------------------------------------------
# High priority programs
LN_BOCO:
	$(CC) $(CFLAGS) -o LN_BOCO src/LN_BOCO.cpp $(LIBRARIES)

LN_MP2RAGE_DNOISE:
	$(CC) $(CFLAGS) -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp $(LIBRARIES)

LN_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp $(LIBRARIES)

LN2_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp $(LIBRARIES)

# -----------------------------------------------------------------------------
# Low priority programs
LN_CORREL2FILES:
	$(CC) $(CFLAGS) -o LN_CORREL2FILES src/LN_CORREL2FILES.cpp $(LIBRARIES)

LN_DIRECT_SMOOTH:
	$(CC) $(CFLAGS) -o LN_DIRECT_SMOOTH src/LN_DIRECT_SMOOTH.cpp $(LIBRARIES)

LN_GRADSMOOTH:
	$(CC) $(CFLAGS) -o LN_GRADSMOOTH src/LN_GRADSMOOTH.cpp $(LIBRARIES)

LN_ZOOM:
	$(CC) $(CFLAGS) -o LN_ZOOM src/LN_ZOOM.cpp $(LIBRARIES)

LN_FLOAT_ME:
	$(CC) $(CFLAGS) -o LN_FLOAT_ME src/LN_FLOAT_ME.cpp $(LIBRARIES)

LN_SHORT_ME:
	$(CC) $(CFLAGS) -o LN_SHORT_ME src/LN_SHORT_ME.cpp $(LIBRARIES)

LN_INT_ME:
	$(CC) $(CFLAGS) -o LN_INT_ME src/LN_INT_ME.cpp $(LIBRARIES)

LN_3DCOLUMNS:
	$(CC) $(CFLAGS) -o LN_3DCOLUMNS src/LN_3DCOLUMNS.cpp $(LIBRARIES)

LN_COLUMNAR_DIST:
	$(CC) $(CFLAGS) -o LN_COLUMNAR_DIST src/LN_COLUMNAR_DIST.cpp $(LIBRARIES)

LN_EXTREMETR:
	$(CC) $(CFLAGS) -o LN_EXTREMETR src/LN_EXTREMETR.cpp $(LIBRARIES)

LN_GFACTOR:
	$(CC) $(CFLAGS) -o LN_GFACTOR src/LN_GFACTOR.cpp $(LIBRARIES)

LN_GROW_LAYERS:
	$(CC) $(CFLAGS) -o LN_GROW_LAYERS src/LN_GROW_LAYERS.cpp $(LIBRARIES)

LN_IMAGIRO:
	$(CC) $(CFLAGS) -o LN_IMAGIRO src/LN_IMAGIRO.cpp $(LIBRARIES)

LN_INTPRO:
	$(CC) $(CFLAGS) -o LN_INTPRO src/LN_INTPRO.cpp $(LIBRARIES)

LN_LEAKY_LAYERS:
	$(CC) $(CFLAGS) -o LN_LEAKY_LAYERS src/LN_LEAKY_LAYERS.cpp $(LIBRARIES)

LN_NOISEME:
	$(CC) $(CFLAGS) -o LN_NOISEME src/LN_NOISEME.cpp $(LIBRARIES)

LN_RAGRUG:
	$(CC) $(CFLAGS) -o LN_RAGRUG src/LN_RAGRUG.cpp $(LIBRARIES)

LN_SKEW:
	$(CC) $(CFLAGS) -o LN_SKEW src/LN_SKEW.cpp $(LIBRARIES)

LN_TEMPSMOOTH:
	$(CC) $(CFLAGS) -o LN_TEMPSMOOTH src/LN_TEMPSMOOTH.cpp $(LIBRARIES)

LN_TRIAL:
	$(CC) $(CFLAGS) -o LN_TRIAL src/LN_TRIAL.cpp $(LIBRARIES)

LN_LOITUMA:
	$(CC) $(CFLAGS) -o LN_LOITUMA src/LN_LOITUMA.cpp $(LIBRARIES)

LN_NOISE_KERNEL:
	$(CC) $(CFLAGS) -o LN_NOISE_KERNEL src/LN_NOISE_KERNEL.cpp $(LIBRARIES)

LN2_DEVEIN:
	$(CC) $(CFLAGS) -o LN2_DEVEIN src/LN2_DEVEIN.cpp $(LIBRARIES)

LN_PHYSIO_PARS:
	$(CC) $(CFLAGS) -o LN_PHYSIO_PARS src/LN_PHYSIO_PARS.cpp $(LIBRARIES)

clean:
	$(RM) obj/*.o $(LAYNII)
