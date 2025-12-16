Original Pinball Compilation

This project is planned to contain multiple pinball tables. As of this writing, only one is available. It is incomplete, but playable.

Build (macOS)

Prereqs:
- Xcode Command Line Tools (`xcode-select --install`)
- Homebrew deps: `brew install sdl2 sdl2_image sdl2_mixer pkg-config`

Option A (script):
- `bash tools/build-macos.sh`
- Run: `bash tools/run.sh`

Option B (CMake):
- `cmake -S . -B build && cmake --build build`
- Run: `./build/pinball`

Option C (Make):
- `make`
- Run: `make run`

Notes:
- The game loads assets via relative paths; the executable will try to auto-select the right working directory at startup.
- You can override the data directory with `OPC_DATA_DIR=/path/to/repo/root ./build/pinball`.

Keyboard controls are as follows:

Start - C

Coin - V

Left Flipper - F

Right Flipper - J

Plunger - M

Left Nudge - D

Right Nudge - K

Front Nudge - Space

Camera Pan - Numpad 8462

Camera Zoom - +/-

Camera Home - Numpad 5

Table rules (Magical Robot):
- Complete A-B-C to light captive ball for double bonus
- Complete drop targets to light bumper (for 100) and outlane (for 5000) on corresponding side
- Complete both drop target banks to light eject hole for extra ball
- Drop targets and eject hole advance bonus by 1000 (max 15000)
- Tilt forfeits ball in play
