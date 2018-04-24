
# note the TARFILE_NAME embeds the release version number
TARFILE_NAME	= nifti2clib-0.0.1

USEZLIB         = -DHAVE_ZLIB

## Compiler  defines
CC		= g++ 
IFLAGS          = -I. 
CFLAGS          = -Wall -pedantic $(USEZLIB) $(IFLAGS)

LDLIBS=     = -gsl
LLIBS 		= -lz -lm -lgsl -lgslcblas 

MISC_OBJS	= nifticdf.o znzlib.o
OBJS	   	= nifti2_io.o $(MISC_OBJS)

TOOLS 	   	= nifti_tool nifti1_tool nifti2_tool

# include your future program below: EXAMPLES   	= My_nii_read    My_future_program_name
EXAMPLES   	= My_nii_read LN_SNR_int_corr LN_FAsim LN_NOISEME LN_GROW_LAYERS LN_3DGROW_LAYERS LN_DEBUGGING LN_GFACTOR LN_LEAKY_LAYERS LN_LAYER_SMOOTH LN_3DCOLUMNS LN_SHORT_ME LN_FAsim LN_FIX_RIM LN_FLOAT_ME LN_IMAGIRO
# main targets (primary is nifti_tool, for now)
nifti_tool: nifti_tool.o nifti_tool.h nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

#all: $(TOOLS) $(EXAMPLES)
all:  $(EXAMPLES)


nifti2objs: $(OBJS)

LN_SNR_int_corr: LN_SNR_int_corr.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

My_nii_read: My_nii_read.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)
	
LN_FAsim: LN_FAsim.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_NOISEME: LN_NOISEME.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_GROW_LAYERS: LN_GROW_LAYERS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)
	
LN_3DGROW_LAYERS: LN_3DGROW_LAYERS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)
	
LN_DEBUGGING: LN_DEBUGGING.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)
	
LN_GFACTOR: LN_GFACTOR.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_LEAKY_LAYERS: LN_LEAKY_LAYERS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_LAYER_SMOOTH: LN_LAYER_SMOOTH.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_3DCOLUMNS: LN_3DCOLUMNS.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_SHORT_ME: LN_SHORT_ME.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)
	
LN_FLOAT_ME: LN_FLOAT_ME.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_FIX_RIM: LN_FIX_RIM.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

LN_IMAGIRO: LN_IMAGIRO.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $< $(OBJS) $(LLIBS)

clean:
	$(RM) *.o $(TOOLS) $(EXAMPLES)

