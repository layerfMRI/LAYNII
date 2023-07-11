# LayNii makefile

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
					LN_INFO \
					LN_CONLAY \
					LN2_DIRECTIONALITY_BIN \


LAYNII2		= 	LN2_LAYERS \
				LN2_COLUMNS \
				LN2_CONNECTED_CLUSTERS \
				LN2_MULTILATERATE \
				LN2_PATCH_FLATTEN \
				LN2_PATCH_UNFLATTEN \
				LN2_PATCH_FLATTEN_2D \
				LN2_CHOLMO \
				LN2_PROFILE \
				LN2_LAYERDIMENSION \
				LN2_MASK \
				LN2_BORDERIZE \
				LN2_GEODISTANCE \
				LN2_IFPOINTS \
				LN2_DEVEIN \
				LN2_RIMIFY \
				LN2_VORONOI \
				LN2_ZERO_CROSSING \
				LN2_UVD_FILTER \
				LN2_HEXBIN \
				LN2_NEIGHBORS \
				LN2_GRAMAG \
				LN2_RIM_POLISH \


LAYNII 	= $(LAYNII2) $(HIGH_PRIORITY) $(LOW_PRIORITY)

# =============================================================================
all : $(LAYNII)

.PHONY: all $(LAYNII2) $(HIGH_PRIORITY) $(LOW_PRIORITY)

# =============================================================================
# LAYNII programs
LN2_LAYERS:
	$(CC) $(CFLAGS) -o LN2_LAYERS src/LN2_LAYERS.cpp $(LIBRARIES) $(LFLAGS)

LN2_LLOYD:
	$(CC) $(CFLAGS) -o obj/LN2_LLOYD.o src/LN2_LLOYD.cpp $(LIBRARIES) $(LFLAGS)

LN_BOCO:
	$(CC) $(CFLAGS) -o LN_BOCO src/LN_BOCO.cpp $(LIBRARIES) $(LFLAGS)

LN_MP2RAGE_DNOISE:
	$(CC) $(CFLAGS) -o LN_MP2RAGE_DNOISE src/LN_MP2RAGE_DNOISE.cpp $(LIBRARIES) $(LFLAGS)

LN_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN_LAYER_SMOOTH src/LN_LAYER_SMOOTH.cpp $(LIBRARIES) $(LFLAGS)

LN2_LAYER_SMOOTH:
	$(CC) $(CFLAGS) -o LN2_LAYER_SMOOTH src/LN2_LAYER_SMOOTH.cpp $(LIBRARIES) $(LFLAGS)

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

LN2_DIRECTIONALITY_BIN:
	$(CC) $(CFLAGS) -o LN2_DIRECTIONALITY_BIN src/LN2_DIRECTIONALITY_BIN.cpp $(LIBRARIES) $(LFLAGS)

LN2_DEVEIN:
	$(CC) $(CFLAGS) -o LN2_DEVEIN src/LN2_DEVEIN.cpp $(LIBRARIES) $(LFLAGS)

LN_PHYSIO_PARS:
	$(CC) $(CFLAGS) -o LN_PHYSIO_PARS src/LN_PHYSIO_PARS.cpp $(LIBRARIES) $(LFLAGS)

LN2_RIMIFY:
	$(CC) $(CFLAGS) -o LN2_RIMIFY src/LN2_RIMIFY.cpp $(LIBRARIES) $(LFLAGS)

LN_INFO:
	$(CC) $(CFLAGS) -o LN_INFO src/LN_INFO.cpp $(LIBRARIES) $(LFLAGS)

LN_CONLAY:
	$(CC) $(CFLAGS) -o LN_CONLAY src/LN_CONLAY.cpp $(LIBRARIES) $(LFLAGS)

LN2_COLUMNS:
	$(CC) $(CFLAGS) -o LN2_COLUMNS src/LN2_COLUMNS.cpp $(LIBRARIES) $(LFLAGS)

LN2_CONNECTED_CLUSTERS:
	$(CC) $(CFLAGS) -o LN2_CONNECTED_CLUSTERS src/LN2_CONNECTED_CLUSTERS.cpp $(LIBRARIES) $(LFLAGS)

LN2_MULTILATERATE:
	$(CC) $(CFLAGS) -o LN2_MULTILATERATE src/LN2_MULTILATERATE.cpp $(LIBRARIES) $(LFLAGS)

LN2_PATCH_FLATTEN:
	$(CC) $(CFLAGS) -o LN2_PATCH_FLATTEN src/LN2_PATCH_FLATTEN.cpp $(LIBRARIES) $(LFLAGS)

LN2_PATCH_FLATTEN_2D:
	$(CC) $(CFLAGS) -o LN2_PATCH_FLATTEN_2D src/LN2_PATCH_FLATTEN_2D.cpp $(LIBRARIES) $(LFLAGS)

LN2_PATCH_UNFLATTEN:
	$(CC) $(CFLAGS) -o LN2_PATCH_UNFLATTEN src/LN2_PATCH_UNFLATTEN.cpp $(LIBRARIES) $(LFLAGS)

LN2_CHOLMO:
	$(CC) $(CFLAGS) -o LN2_CHOLMO src/LN2_CHOLMO.cpp $(LIBRARIES) $(LFLAGS)

LN2_PROFILE:
	$(CC) $(CFLAGS) -o LN2_PROFILE src/LN2_PROFILE.cpp $(LIBRARIES) $(LFLAGS)

LN2_LAYERDIMENSION:
	$(CC) $(CFLAGS) -o LN2_LAYERDIMENSION src/LN2_LAYERDIMENSION.cpp $(LIBRARIES) $(LFLAGS)

LN2_MASK:
	$(CC) $(CFLAGS) -o LN2_MASK src/LN2_MASK.cpp $(LIBRARIES) $(LFLAGS)

LN2_BORDERIZE:
	$(CC) $(CFLAGS) -o LN2_BORDERIZE src/LN2_BORDERIZE.cpp $(LIBRARIES) $(LFLAGS)

LN2_GEODISTANCE:
	$(CC) $(CFLAGS) -o LN2_GEODISTANCE src/LN2_GEODISTANCE.cpp $(LIBRARIES) $(LFLAGS)

LN2_IFPOINTS:
	$(CC) $(CFLAGS) -o LN2_IFPOINTS src/LN2_IFPOINTS.cpp $(LIBRARIES) $(LFLAGS)

LN2_VORONOI:
	$(CC) $(CFLAGS) -o LN2_VORONOI src/LN2_VORONOI.cpp $(LIBRARIES) $(LFLAGS)

LN2_ZERO_CROSSING:
	$(CC) $(CFLAGS) -o LN2_ZERO_CROSSING src/LN2_ZERO_CROSSING.cpp $(LIBRARIES) $(LFLAGS)

LN2_HEXBIN:
	$(CC) $(CFLAGS) -o LN2_HEXBIN src/LN2_HEXBIN.cpp $(LIBRARIES) $(LFLAGS)

LN2_UVD_FILTER:
	$(CC) $(CFLAGS) -o LN2_UVD_FILTER src/LN2_UVD_FILTER.cpp $(LIBRARIES) $(LFLAGS)

# =============================================================================
# Work in progress programs
LN2_GRAMAG:
	$(CC) $(CFLAGS) -o LN2_GRAMAG src/LN2_GRAMAG.cpp $(LIBRARIES) $(LFLAGS)

LN2_UVD_LSTSQR:
	$(CC) $(CFLAGS) -o LN2_UVD_LSTSQR src/LN2_UVD_LSTSQR.cpp $(LIBRARIES) $(LFLAGS)

LN2_PEAK_DETECT:
	$(CC) $(CFLAGS) -o LN2_PEAK_DETECT src/LN2_PEAK_DETECT.cpp $(LIBRARIES) $(LFLAGS)

LN2_NEIGHBORS:
	$(CC) $(CFLAGS) -o LN2_NEIGHBORS src/LN2_NEIGHBORS.cpp $(LIBRARIES) $(LFLAGS)

LN2_WINDOWED_COUNTER_2D:
	$(CC) $(CFLAGS) -o LN2_WINDOWED_COUNTER_2D src/LN2_WINDOWED_COUNTER_2D.cpp $(LIBRARIES) $(LFLAGS)

LN2_RIM_POLISH:
	$(CC) $(CFLAGS) -o LN2_RIM_POLISH src/LN2_RIM_POLISH.cpp $(LIBRARIES) $(LFLAGS)

# =============================================================================

clean:
	$(RM) obj/*.o $(LAYNII)

tests:
	cd test_data && bash ./tests.sh

# =============================================================================
# Maintenance
# =============================================================================
bump_version:
	VERSION=`cat dep/laynii_lib.cpp | grep LayNii | cut -c 21- | rev | cut -c 20- | rev` && \
	echo "\nbumping to version: $$VERSION \n" && \
	sed -i "s/^version:.*/version: $$VERSION/g" CITATION.cff
	make Dockerfile


# =============================================================================
# Docker related content
# =============================================================================
.PHONY: Dockerfile docker_build

Dockerfile:
	VERSION=`cat CITATION.cff | grep ^version | cut -c 10-` && \
	docker run --rm repronim/neurodocker:0.7.0 generate docker \
		--base debian:stretch-slim \
		--label version=$$VERSION \
		--pkg-manager apt \
		--install "build-essential libz-dev" \
		--run "mkdir /input /output" \
		--run "chmod -R 777 /output" \
		--user laynii \
		--copy dep /home/laynii/dep \
		--copy src /home/laynii/src \
		--copy Makefile /home/laynii/ \
		--workdir /home/laynii/ \
		--run "make all" \
		--add-to-entrypoint 'PATH=/home/laynii/:$$PATH' > Dockerfile

docker_build:
	VERSION=`cat CITATION.cff | grep ^version | cut -c 10-` && \
	docker build . -t laynii:$$VERSION -t laynii:latest
