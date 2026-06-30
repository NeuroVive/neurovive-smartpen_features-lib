# =============================================================================
# NeuroVive Smart Pen — Build Configuration
# =============================================================================
# Targets:
#   smartpen_test       — Standalone test harness (main.cpp + lib)
#   libsmartpen.so      — Shared library for Dart FFI consumption
#   libsmartpen.a       — Static library (embedded deployment)
#   clean               — Remove all build artifacts
#
# Toolchain requirements: any C++17 compiler (g++ ≥ 8, clang++ ≥ 7).
# =============================================================================

CXX          ?= g++
CXXSTD       ?= -std=c++17
CXXWARN      ?= -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
CXXOPT       ?= -O2 -ffast-math -fno-exceptions -fno-rtti
CXXINC       := -I.

# Emit position-independent code so the same objects feed both the static
# archive and the FFI shared library.
CXXFLAGS     := $(CXXSTD) $(CXXWARN) $(CXXOPT) $(CXXINC) -fPIC

AR           ?= ar
ARFLAGS      := rcs

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
LIB_OBJS     := smartPen_Features.o
TEST_BIN     := smartpen_test
SHARED_LIB   := libsmartpen.so
STATIC_LIB   := libsmartpen.a

.PHONY: all clean test

all: $(TEST_BIN) $(SHARED_LIB) $(STATIC_LIB)

# Standalone test harness
$(TEST_BIN): main.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

# Dart FFI shared library
$(SHARED_LIB): $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ -lm

# Embedded static archive
$(STATIC_LIB): $(LIB_OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Object rules
%.o: %.cpp smartPen_Features.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(LIB_OBJS) main.o $(TEST_BIN) $(SHARED_LIB) $(STATIC_LIB)
