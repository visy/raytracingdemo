g++ -Os -fdata-sections -ffunction-sections -fipa-pta demo.cpp pixie.cpp pixie_win.cpp -D DEMO_FULLSCREEN=true -Wl,--gc-sections -Wl,-O1 -Wl,--as-needed -Wl,--strip-all  -lgdi32 -o demo.exe
upx --best demo.exe
demo.exe