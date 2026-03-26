# Notes for LayNii IDA (Integrated Development/Discovery Application)

- Use `make clean` then `make LayNii_IDA` inside `src` folder.
- `make clean && make LayNii_IDA` for quick clean and compile.
- It is strongly recommended to use `g++` compiler (e.g. `brew install gcc`). This boost optimization results significantly compared to `c++` when using `-O3`.

## Windows compilation notes
- I use msys64->mingw64.exe as terminal in windows.
- Here I pull LayNii using git
- Then compile using `make LayNii_IDA`
- I add the dlls if it complains into the exe folder.
- A few useful command I used: `pacman -S mingw-w64_x86_64-pkgconf`, `pacman -S mingw-w64_x86_64-SDL2`, `pkg-config --cflags sdl2`
