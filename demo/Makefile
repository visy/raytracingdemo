CUDA_PATH     ?= /usr/local/cuda
HOST_COMPILER  = g++
NVCC           = $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

# select one of these for Debug vs. Release
#NVCC_DBG       = -g -G
NVCC_DBG       =

NVCCFLAGS      = $(NVCC_DBG) -m64 -lwsock32 -lws2_32 -fpermissive -D DEMO_FULLSCREEN=false -lgdi32 -fopenmp -lbass -O3 -march=skylake -mno-vzeroupper -s
GENCODE_FLAGS  = -gencode arch=compute_60,code=sm_60

SRCS = demo.cu pixie.cpp pixie_win.cpp device.c track.c  
INCS = vec3.h ray.h hitable.h hitable_list.h sphere.h camera.h material.h

demo.exe: demo.o
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o demo.exe demo.o

cudart.o: $(SRCS) $(INCS)
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o demo.o -c demo.cu


clean:
	rm -f demo.exe demo.o
