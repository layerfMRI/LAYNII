
# Note the TARFILE_NAME embeds the release version number
TARFILE_NAME	= nifti2clib-0.0.1

USEZLIB         = -DHAVE_ZLIB

## Compiler  defines
CC				= g++
IFLAGS          = -I.
CFLAGS          = -Wall -pedantic $(USEZLIB) $(IFLAGS)

##
OBJS	   	= obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o
TOOLS 	   	= nifti_tool nifti1_tool nifti2_tool

LAYNII   =	LN_NOISEME \
			LN_GROW_LAYERS \
			LN_DEBUGGING \
			LN_GFACTOR \
			LN_LEAKY_LAYERS \
			LN_LAYER_SMOOTH \
			LN_3DCOLUMNS \
			LN_FLOAT_ME \
			LN_IMAGIRO \
			LN_DIRECT_SMOOTH \
			LN_RAGRUG \
			LN_CORREL2FILES \
			LN_EXTREMETR \
			LN_BOCO \
			LN_TRIAL \
			LN_ZOOM \
			LN_COLUMNAR_DIST \
			LN_GRADSMOOTH \
			LN_SKEW \
			LN_INTPRO \
			LN_TEMPSMOOTH \
			LN_MP2RAGE_DNOISE \
			LN_PHYSIO_PARS

# main targets (primary is nifti_tool, for now)
nifti_tool: obj/nifti_tool.o nifti_tool.h nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)

all:  $(LAYNII)

nifti2objs: $(OBJS)


LN_NOISEME: obj/LN_NOISEME.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_GROW_LAYERS: obj/LN_GROW_LAYERS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_GFACTOR: obj/LN_GFACTOR.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_LEAKY_LAYERS: obj/LN_LEAKY_LAYERS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_LAYER_SMOOTH: obj/LN_LAYER_SMOOTH.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_3DCOLUMNS: obj/LN_3DCOLUMNS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_FLOAT_ME: obj/LN_FLOAT_ME.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_IMAGIRO: obj/LN_IMAGIRO.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_DIRECT_SMOOTH: obj/LN_DIRECT_SMOOTH.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_RAGRUG: obj/LN_RAGRUG.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_CORREL2FILES: obj/LN_CORREL2FILES.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_EXTREMETR: obj/LN_EXTREMETR.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_BOCO: obj/LN_BOCO.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_TRIAL: obj/LN_TRIAL.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_ZOOM: obj/LN_ZOOM.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_COLUMNAR_DIST: obj/LN_COLUMNAR_DIST.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_GRADSMOOTH: obj/LN_GRADSMOOTH.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_SKEW: obj/LN_SKEW.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_INTPRO: obj/LN_INTPRO.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_TEMPSMOOTH: obj/LN_TEMPSMOOTH.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_MP2RAGE_DNOISE: obj/LN_MP2RAGE_DNOISE.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
LN_PHYSIO_PARS: obj/LN_PHYSIO_PARS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)

clean:
	$(RM) *.o $(TOOLS) $(LAYNII)
