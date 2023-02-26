g++ demo.cpp pixie.cpp pixie_win.cpp device.c track.c -lwsock32 -lws2_32 -fpermissive -D DEMO_FULLSCREEN=true -lgdi32 -fopenmp -lbass -O3 -s -o demo.exe
demo.exe
