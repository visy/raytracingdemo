g++ demo.cpp pixie.cpp pixie_win.cpp device.c track.c -lwsock32 -lws2_32 -fpermissive -D DEMO_FULLSCREEN=false -lgdi32 -fopenmp -lbass -O3 -march=skylake -mno-vzeroupper -mwindows -s -o demo.exe
demo.exe
