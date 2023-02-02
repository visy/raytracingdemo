g++ -s -Os -fdata-sections -ffunction-sections -fipa-pta demo.cpp pixie.cpp pixie_win.cpp -D DEMO_FULLSCREEN=true -Os -s -fopenmp -lbass -lgdi32 -o demo.exe
strip -g demo.exe
demo.exe
