# LAYNII makefile

CC		= g++
IFLAGS	= -I.
USEZLIB	= -DHAVE_ZLIB
CFLAGS	= -std=c++11 -Wall -pedantic $(USEZLIB) $(IFLAGS)

# =============================================================================
LIBRARIES		=	obj/nifti2_io.o \
 					obj/nifticdf.o \
 					obj/znzlib.o \
 					obj/laynii_lib.o \

HIGH_PRIORITY	= 	LN_BOCO \
					LN_MP2RAGE_DNOISE \
					LN_LAYER_SMOOTH \

LOW_PRIORITY	=	LN_3DCOLUMNS \
					LN_COLUMNAR_DIST \
					LN_CORREL2FILES \
					LN_DIRECT_SMOOTH \
					LN_GRADSMOOTH \
					LN_ZOOM \
					LN_FLOAT_ME \
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

LAYNII2	= LN2_LAYERS

LAYNII 	= $(HIGH_PRIORITY) $(LOW_PRIORITY) $(LAYNII2)

# =============================================================================
dependencies :
	$(CC) -c -std=c++11 -o obj/nifti2_io.o dep/nifti2_io.cpp
	$(CC) -c -std=c++11 -o obj/nifticdf.o dep/nifticdf.cpp
	$(CC) -c -std=c++11 -o obj/znzlib.o dep/znzlib.cpp
	$(CC) -c -std=c++11 -o obj/laynii_lib.o dep/laynii_lib.cpp

all : dependencies $(LAYNII)

# =============================================================================
# LAYNII v2.0.0 programs
LN2_LAYERS:
	$(CC) -c -o -std=c++11 obj/LN2_LAYERS.o src/LN2_LAYERS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN2_LAYERS.o $(LIBRARIES)

# -----------------------------------------------------------------------------
# High priority programs
LN_BOCO:
	$(CC) -c -o -std=c++11 obj/LN_BOCO.o src/LN_BOCO.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_BOCO.o $(LIBRARIES)

LN_MP2RAGE_DNOISE:
	$(CC) -c -o -std=c++11 obj/LN_MP2RAGE_DNOISE.o src/LN_MP2RAGE_DNOISE.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_MP2RAGE_DNOISE.o $(LIBRARIES)

LN_LAYER_SMOOTH:
	$(CC) -c -o -std=c++11 obj/LN_LAYER_SMOOTH.o src/LN_LAYER_SMOOTH.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_LAYER_SMOOTH.o $(LIBRARIES)

# -----------------------------------------------------------------------------
# Low priority programs
LN_CORREL2FILES:
	$(CC) -c -o -std=c++11 obj/LN_BOCO.o src/LN_BOCO.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_BOCO.o $(LIBRARIES)

LN_DIRECT_SMOOTH:
	$(CC) -c -o -std=c++11 obj/LN_DIRECT_SMOOTH.o src/LN_DIRECT_SMOOTH.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_DIRECT_SMOOTH.o $(LIBRARIES)

LN_GRADSMOOTH:
	$(CC) -c -o -std=c++11 obj/LN_GRADSMOOTH.o src/LN_GRADSMOOTH.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_GRADSMOOTH.o $(LIBRARIES)

LN_ZOOM:
	$(CC) -c -o -std=c++11 obj/LN_ZOOM.o src/LN_ZOOM.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_ZOOM.o $(LIBRARIES)

LN_FLOAT_ME:
	$(CC) -c -o -std=c++11 obj/LN_FLOAT_ME.o src/LN_FLOAT_ME.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_FLOAT_ME.o $(LIBRARIES)

LN_3DCOLUMNS:
	$(CC) -c -o -std=c++11 obj/LN_3DCOLUMNS.o src/LN_3DCOLUMNS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_3DCOLUMNS.o $(LIBRARIES)

LN_COLUMNAR_DIST:
	$(CC) -c -o -std=c++11 obj/LN_COLUMNAR_DIST.o src/LN_COLUMNAR_DIST.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_COLUMNAR_DIST.o $(LIBRARIES)

LN_EXTREMETR:
	$(CC) -c -o -std=c++11 obj/LN_EXTREMETR.o src/LN_EXTREMETR.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_EXTREMETR.o $(LIBRARIES)

LN_GFACTOR:
	$(CC) -c -o -std=c++11 obj/LN_GFACTOR.o src/LN_GFACTOR.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_GFACTOR.o $(LIBRARIES)

LN_GROW_LAYERS:
	$(CC) -c -o -std=c++11 obj/LN_GROW_LAYERS.o src/LN_GROW_LAYERS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_GROW_LAYERS.o $(LIBRARIES)

LN_IMAGIRO:
	$(CC) -c -o -std=c++11 obj/LN_IMAGIRO.o src/LN_IMAGIRO.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_IMAGIRO.o $(LIBRARIES)

LN_INTPRO:
	$(CC) -c -o -std=c++11 obj/LN_INTPRO.o src/LN_INTPRO.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_INTPRO.o $(LIBRARIES)

LN_LEAKY_LAYERS:
	$(CC) -c -o -std=c++11 obj/LN_LEAKY_LAYERS.o src/LN_LEAKY_LAYERS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_LEAKY_LAYERS.o $(LIBRARIES)

LN_NOISEME:
	$(CC) -c -o -std=c++11 obj/LN_NOISEME.o src/LN_NOISEME.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_NOISEME.o $(LIBRARIES)

LN_RAGRUG:
	$(CC) -c -o -std=c++11 obj/LN_RAGRUG.o src/LN_RAGRUG.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_RAGRUG.o $(LIBRARIES)

LN_SKEW:
	$(CC) -c -o -std=c++11 obj/LN_SKEW.o src/LN_SKEW.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_SKEW.o $(LIBRARIES)

LN_TEMPSMOOTH:
	$(CC) -c -o -std=c++11 obj/LN_TEMPSMOOTH.o src/LN_TEMPSMOOTH.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_TEMPSMOOTH.o $(LIBRARIES)

LN_TRIAL:
	$(CC) -c -o -std=c++11 obj/LN_TRIAL.o src/LN_TRIAL.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_TRIAL.o $(LIBRARIES)

LN_PHYSIO_PARS:
	$(CC) -c -o -std=c++11 obj/LN_PHYSIO_PARS.o src/LN_PHYSIO_PARS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_PHYSIO_PARS.o $(LIBRARIES)

clean:
	$(RM) obj/*.o $(LAYNII)
