# LAYNII makefile

CC		= c++
CFLAGS	= -std=c++11 -DHAVE_ZLIB 
LFLAGS	= -lm -lz
# CFLAGS	= -std=c++11 -pedantic -DHAVE_ZLIB -lm -lz

# =============================================================================
LIBRARIES		=	dep/nifti2_io.cpp \
					dep/znzlib.cpp \
					dep/laynii_lib.cpp \
					-I./dep \

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
	$(CC) $(CFLAGS) -o LN2_LAYERS src/LN2_LAYERS.cpp $(LIBRARIES) $(LFLAGS)

LN2_LLOYD:
	$(CC) $(CFLAGS) -o obj/LN2_LLOYD.o src/LN2_LLOYD.cpp $(LIBRARIES) $(LFLAGS)

# -----------------------------------------------------------------------------
# High priority programs
LN_BOCO:
	$(CC) $(CFLAGS) -o LN_BOCO src/LN_BOCO.cpp $(LIBRARIES) $(LFLAGS)

LN_MP2RAGE_DNOISE:
	$(CC) $(CFLAGS) -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp $(LIBRARIES) $(LFLAGS)

LN_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp $(LIBRARIES) $(LFLAGS)

LN2_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp $(LIBRARIES) $(LFLAGS)

# -----------------------------------------------------------------------------
# Low priority programs
LN_CORREL2FILES:
	$(CC) $(CFLAGS) -o LN_CORREL2FILES src/LN_CORREL2FILES.cpp $(LIBRARIES) $(LFLAGS)

LN_DIRECT_SMOOTH:
	$(CC) $(CFLAGS) -o LN_DIRECT_SMOOTH src/LN_DIRECT_SMOOTH.cpp $(LIBRARIES) $(LFLAGS)

LN_GRADSMOOTH:
	$(CC) $(CFLAGS) -o LN_GRADSMOOTH src/LN_GRADSMOOTH.cpp $(LIBRARIES) $(LFLAGS)

LN_ZOOM:
	$(CC) $(CFLAGS) -o LN_ZOOM src/LN_ZOOM.cpp $(LIBRARIES) $(LFLAGS)

LN_FLOAT_ME:
	$(CC) $(CFLAGS) -o LN_FLOAT_ME src/LN_FLOAT_ME.cpp $(LIBRARIES) $(LFLAGS)

LN_SHORT_ME:
	$(CC) $(CFLAGS) -o LN_SHORT_ME src/LN_SHORT_ME.cpp $(LIBRARIES) $(LFLAGS)

LN_INT_ME:
	$(CC) $(CFLAGS) -o LN_INT_ME src/LN_INT_ME.cpp $(LIBRARIES) $(LFLAGS)

LN_3DCOLUMNS:
	$(CC) $(CFLAGS) -o LN_3DCOLUMNS src/LN_3DCOLUMNS.cpp $(LIBRARIES) $(LFLAGS)

LN_COLUMNAR_DIST:
	$(CC) $(CFLAGS) -o LN_COLUMNAR_DIST src/LN_COLUMNAR_DIST.cpp $(LIBRARIES) $(LFLAGS)

LN_EXTREMETR:
	$(CC) $(CFLAGS) -o LN_EXTREMETR src/LN_EXTREMETR.cpp $(LIBRARIES) $(LFLAGS)

LN_GFACTOR:
	$(CC) $(CFLAGS) -o LN_GFACTOR src/LN_GFACTOR.cpp $(LIBRARIES) $(LFLAGS)

LN_GROW_LAYERS:
	$(CC) $(CFLAGS) -o LN_GROW_LAYERS src/LN_GROW_LAYERS.cpp $(LIBRARIES) $(LFLAGS)

LN_IMAGIRO:
	$(CC) $(CFLAGS) -o LN_IMAGIRO src/LN_IMAGIRO.cpp $(LIBRARIES) $(LFLAGS)

LN_INTPRO:
	$(CC) $(CFLAGS) -o LN_INTPRO src/LN_INTPRO.cpp $(LIBRARIES) $(LFLAGS)

LN_LEAKY_LAYERS:
	$(CC) $(CFLAGS) -o LN_LEAKY_LAYERS src/LN_LEAKY_LAYERS.cpp $(LIBRARIES) $(LFLAGS)

LN_NOISEME:
	$(CC) $(CFLAGS) -o LN_NOISEME src/LN_NOISEME.cpp $(LIBRARIES) $(LFLAGS)

LN_RAGRUG:
	$(CC) $(CFLAGS) -o LN_RAGRUG src/LN_RAGRUG.cpp $(LIBRARIES) $(LFLAGS)

LN_SKEW:
	$(CC) $(CFLAGS) -o LN_SKEW src/LN_SKEW.cpp $(LIBRARIES) $(LFLAGS)

LN_TEMPSMOOTH:
	$(CC) $(CFLAGS) -o LN_TEMPSMOOTH src/LN_TEMPSMOOTH.cpp $(LIBRARIES) $(LFLAGS)

LN_TRIAL:
	$(CC) $(CFLAGS) -o LN_TRIAL src/LN_TRIAL.cpp $(LIBRARIES) $(LFLAGS)

LN_LOITUMA:
	$(CC) $(CFLAGS) -o LN_LOITUMA src/LN_LOITUMA.cpp $(LIBRARIES) $(LFLAGS)

LN_NOISE_KERNEL:
	$(CC) $(CFLAGS) -o LN_NOISE_KERNEL src/LN_NOISE_KERNEL.cpp $(LIBRARIES) $(LFLAGS)

LN2_DEVEIN:
	$(CC) $(CFLAGS) -o LN2_DEVEIN src/LN2_DEVEIN.cpp $(LIBRARIES) $(LFLAGS)

LN_PHYSIO_PARS:
	$(CC) $(CFLAGS) -o LN_PHYSIO_PARS src/LN_PHYSIO_PARS.cpp $(LIBRARIES) $(LFLAGS)

clean:
	$(RM) obj/*.o $(LAYNII)
