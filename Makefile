CXX ?= clang++
PKG_CONFIG ?= pkg-config

SDL_PKGS ?= sdl2 SDL2_image SDL2_mixer
CPPFLAGS += $(shell $(PKG_CONFIG) --cflags $(SDL_PKGS))
LDLIBS += $(shell $(PKG_CONFIG) --libs $(SDL_PKGS))

CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17

BUILD_DIR := build
TARGET := $(BUILD_DIR)/pinball

.PHONY: all clean run

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): engine.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR)

