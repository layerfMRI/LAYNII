# LAYNII makefile

CC		= g++
IFLAGS	= -I.
USEZLIB	= -DHAVE_ZLIB
CFLAGS	= -Wall -pedantic $(USEZLIB) $(IFLAGS)

LIBRARIES	= obj/nifti2_io.o
LIBRARIES	+= obj/nifticdf.o
LIBRARIES	+= obj/znzlib.o
LIBRARIES	+= obj/laynii_lib.o

HIGH_PRIORITY	= LN_BOCO
HIGH_PRIORITY	+= LN_MP2RAGE_DNOISE
HIGH_PRIORITY	+= LN_LAYER_SMOOTH

LAYNII2		= LN2_GROW_LAYERS

LAYNII	=	$(HIGH_PRIORITY) \
			# LN_3DCOLUMNS \
			# LN_COLUMNAR_DIST \
			# LN_CORREL2FILES \
			# LN_DEBUGGING \
			# LN_DIRECT_SMOOTH \
			# LN_EXTREMETR \
			# LN_FLOAT_ME \
			# LN_GFACTOR \
			# LN_GRADSMOOTH \
			# LN_GROW_LAYERS \
			# LN_IMAGIRO \
			# LN_INTPRO \
			# LN_LAYER_SMOOTH \
			# LN_LEAKY_LAYERS \
			# LN_NOISEME \
			# LN_PHYSIO_PARS \
			# LN_RAGRUG \
			# LN_SKEW \
			# LN_TEMPSMOOTH \
			# LN_TRIAL \
			# LN_ZOOM

dependencies :
	$(CC) -c -o obj/nifti2_io.o dep/nifti2_io.cpp
	$(CC) -c -o obj/nifticdf.o dep/nifticdf.cpp
	$(CC) -c -o obj/znzlib.o dep/znzlib.cpp
	$(CC) -c -o obj/laynii_lib.o dep/laynii_lib.cpp

all : dependencies $(LAYNII) $(LAYNII2)

# -----------------------------------------------------------------------------
# LAYNII v2.0.0 programs
LN2_GROW_LAYERS:
	$(CC) -c -o obj/LN2_GROW_LAYERS.o src/LN2_GROW_LAYERS.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN2_GROW_LAYERS.o $(LIBRARIES)

# -----------------------------------------------------------------------------
# High priority programs
LN_BOCO:
	$(CC) -c -o obj/LN_BOCO.o src/LN_BOCO.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_BOCO.o $(LIBRARIES)

LN_MP2RAGE_DNOISE:
	$(CC) -c -o obj/LN_MP2RAGE_DNOISE.o src/LN_MP2RAGE_DNOISE.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_MP2RAGE_DNOISE.o $(LIBRARIES)

LN_LAYER_SMOOTH:
	$(CC) -c -o obj/LN_LAYER_SMOOTH.o src/LN_LAYER_SMOOTH.cpp
	$(CC) -o $@ $(CFLAGS) obj/LN_LAYER_SMOOTH.o $(LIBRARIES)

# -----------------------------------------------------------------------------
# Other programs
# LN_3DCOLUMNS: obj/LN_3DCOLUMNS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_COLUMNAR_DIST: obj/LN_COLUMNAR_DIST.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_CORREL2FILES: obj/LN_CORREL2FILES.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_DIRECT_SMOOTH: obj/LN_DIRECT_SMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_EXTREMETR: obj/LN_EXTREMETR.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_FLOAT_ME: obj/LN_FLOAT_ME.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_GFACTOR: obj/LN_GFACTOR.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_GRADSMOOTH: obj/LN_GRADSMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_GROW_LAYERS: obj/LN_GROW_LAYERS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_IMAGIRO: obj/LN_IMAGIRO.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_INTPRO: obj/LN_INTPRO.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_LEAKY_LAYERS: obj/LN_LEAKY_LAYERS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_NOISEME: obj/LN_NOISEME.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_PHYSIO_PARS: obj/LN_PHYSIO_PARS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_RAGRUG: obj/LN_RAGRUG.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_SKEW: obj/LN_SKEW.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_TEMPSMOOTH: obj/LN_TEMPSMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_TRIAL: obj/LN_TRIAL.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)
# LN_ZOOM: obj/LN_ZOOM.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(LIBRARIES)

clean:
	$(RM) obj/*.o $(LAYNII)
