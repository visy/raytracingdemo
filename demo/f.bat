g++ -s -Os -fdata-sections -ffunction-sections -fipa-pta demo.cpp pixie.cpp pixie_win.cpp -D DEMO_FULLSCREEN=true -Wl,--gc-sections -Wl,-O1 -Wl,--as-needed -Wl,--strip-all -lbass -lgdi32 -o demo.exe
strip -g demo.exe
demo.exe
