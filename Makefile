## Compiler defines
CC		= g++
IFLAGS	= -I.
USEZLIB	= -DHAVE_ZLIB
CFLAGS	= -Wall -pedantic $(USEZLIB) $(IFLAGS)

##
OBJS	= obj/nifti2_io.o obj/nifticdf.o obj/znzlib.o obj/laynii_lib.o
TOOLS	= nifti_tool nifti1_tool nifti2_tool

LAYNII	=	LN2_GROW_LAYERS \
			# LN_3DCOLUMNS \
			# LN_BOCO \
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
			# LN_MP2RAGE_DNOISE \
			# LN_NOISEME \
			# LN_PHYSIO_PARS \
			# LN_RAGRUG \
			# LN_SKEW \
			# LN_TEMPSMOOTH \
			# LN_TRIAL \
			# LN_ZOOM \

# main targets (primary is nifti_tool, for now)
nifti_tool: obj/nifti_tool.o nifti_tool.h objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)

all:  $(OBJS) $(LAYNII)

objs: $(OBJS)

LN2_GROW_LAYERS: obj/LN2_GROW_LAYERS.o objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS)

# LN_3DCOLUMNS: obj/LN_3DCOLUMNS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_BOCO: obj/LN_BOCO.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_COLUMNAR_DIST: obj/LN_COLUMNAR_DIST.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_CORREL2FILES: obj/LN_CORREL2FILES.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_DIRECT_SMOOTH: obj/LN_DIRECT_SMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_EXTREMETR: obj/LN_EXTREMETR.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_FLOAT_ME: obj/LN_FLOAT_ME.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_GFACTOR: obj/LN_GFACTOR.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_GRADSMOOTH: obj/LN_GRADSMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_GROW_LAYERS: obj/LN_GROW_LAYERS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_IMAGIRO: obj/LN_IMAGIRO.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_INTPRO: obj/LN_INTPRO.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_LAYER_SMOOTH: obj/LN_LAYER_SMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_LEAKY_LAYERS: obj/LN_LEAKY_LAYERS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_MP2RAGE_DNOISE: obj/LN_MP2RAGE_DNOISE.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_NOISEME: obj/LN_NOISEME.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_PHYSIO_PARS: obj/LN_PHYSIO_PARS.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_RAGRUG: obj/LN_RAGRUG.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_SKEW: obj/LN_SKEW.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_TEMPSMOOTH: obj/LN_TEMPSMOOTH.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_TRIAL: obj/LN_TRIAL.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)
# LN_ZOOM: obj/LN_ZOOM.o objs
# 	$(CC) -o $@ $(CFLAGS) $< $(OBJS)

clean:
	$(RM) obj/*.o $(TOOLS) $(LAYNII)
