Original Pinball Compilation

This project is planned to contain multiple pinball tables. As of this writing, two are available. One is complete; the other is incomplete, but playable.

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

Left Secondary - G

Right Flipper - J

Right Secondary - H (currently unused)

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

Table rules (Saturn Strikes):
- Spinner lights 5-bank target for Strike or Spare
- Hit lit 5-bank target to advance bowling multiplier and bonus multiplier and score 200k for Strike or 100k for Spare (subject to bowling multiplier)
- Completing 5-bank without hitting lit target only advances bonus multiplier
- Complete 1-2-3 to advance bonus multiplier (lane change on right flipper)
- Left secondary button activates outlane saver magnet when lit
- Standup target lights spinner for 10x score and double bonus on next rip
- Skillshot on standup target scores 10k
- Hit strobing 3-bank target to advance eject hole reward (10k, 30k, 50k, collect bonus) and score 10k
- Complete 3-bank to relight Hook Shot
- Strike with 3x bowling multiplier lights extra ball on outlanes, then special on standup target
- Rollovers, targets, and spinner advance bonus
- Super bonus held at 20k, 30k, and 40k
- Tilt forfeits ball in play
