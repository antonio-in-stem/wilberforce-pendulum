# =============================================================================
#  Wilberforce Pendulum — GNU Make build
#
#  Targets:
#     make            build the solver           -> build/wilberforce(.exe)
#     make run        build and run preset 14
#     make test       build and run the test suite
#     make plot       build, run, and render PNGs with Gnuplot
#     make clean      remove build artefacts
#
#  Works with GNU make + gfortran on Linux/macOS and mingw32-make on Windows.
# =============================================================================

FC      := gfortran
FFLAGS  := -O2 -std=f2008 -Wall -Wextra -fimplicit-none
BUILD   := build
SRCDIR  := src
APPDIR  := app
TESTDIR := test

# Executable suffix (.exe on Windows)
ifeq ($(OS),Windows_NT)
    EXE := .exe
    RM  := del /Q /F
else
    EXE :=
    RM  := rm -f
endif

# Module sources IN DEPENDENCY ORDER (order matters for .mod generation).
MODS := \
    $(SRCDIR)/wilberforce_kinds.f90 \
    $(SRCDIR)/wilberforce_types.f90 \
    $(SRCDIR)/wilberforce_model.f90 \
    $(SRCDIR)/wilberforce_analytic.f90 \
    $(SRCDIR)/wilberforce_energy.f90 \
    $(SRCDIR)/wilberforce_io.f90 \
    $(SRCDIR)/wilberforce_integrator.f90 \
    $(SRCDIR)/wilberforce_presets.f90 \
    $(SRCDIR)/wilberforce_gnuplot.f90

OBJS := $(patsubst $(SRCDIR)/%.f90,$(BUILD)/%.o,$(MODS))
EXEC := $(BUILD)/wilberforce$(EXE)
TEXEC:= $(BUILD)/test_wilberforce$(EXE)

.PHONY: all run test plot clean dirs

all: $(EXEC)

dirs:
	@mkdir -p $(BUILD)

# Compile modules in order. -J places .mod files in $(BUILD).
$(BUILD)/%.o: $(SRCDIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(EXEC): $(OBJS) $(APPDIR)/main.f90
	$(FC) $(FFLAGS) -J$(BUILD) $(OBJS) $(APPDIR)/main.f90 -o $@

$(TEXEC): $(OBJS) $(TESTDIR)/test_wilberforce.f90
	$(FC) $(FFLAGS) -J$(BUILD) $(OBJS) $(TESTDIR)/test_wilberforce.f90 -o $@

run: $(EXEC)
	$(EXEC) 14

plot: $(EXEC)
	$(EXEC) 14 --plot

test: $(TEXEC)
	$(TEXEC)

clean:
	-$(RM) $(BUILD)/*.o $(BUILD)/*.mod $(EXEC) $(TEXEC)

# --- explicit module dependencies (so a change forces the right rebuilds) -----
$(BUILD)/wilberforce_types.o:      $(BUILD)/wilberforce_kinds.o
$(BUILD)/wilberforce_model.o:      $(BUILD)/wilberforce_types.o
$(BUILD)/wilberforce_analytic.o:   $(BUILD)/wilberforce_types.o
$(BUILD)/wilberforce_energy.o:     $(BUILD)/wilberforce_types.o
$(BUILD)/wilberforce_io.o:         $(BUILD)/wilberforce_energy.o
$(BUILD)/wilberforce_integrator.o: $(BUILD)/wilberforce_io.o $(BUILD)/wilberforce_model.o $(BUILD)/wilberforce_analytic.o
$(BUILD)/wilberforce_presets.o:    $(BUILD)/wilberforce_types.o
$(BUILD)/wilberforce_gnuplot.o:    $(BUILD)/wilberforce_io.o
