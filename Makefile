# Compiler settings
CC = gcc
FC = gfortran
CFLAGS = -O3 -march=native -funroll-loops -ffast-math -flto -Wall -g
FFLAGS = -O2 -g
LDFLAGS = -lm

# Directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

SRC_C_DIR = src_c
SRC_F_DIR = src_fortran

# Source files
C_SRC = $(SRC_C_DIR)/main.c $(SRC_C_DIR)/stats.c $(SRC_C_DIR)/data_io.c $(SRC_C_DIR)/mcmc.c
F_SRC = $(SRC_F_DIR)/RandomDistributions.f90 \
        $(SRC_F_DIR)/mod_types.f90 \
        $(SRC_F_DIR)/mod_io.f90 \
        $(SRC_F_DIR)/mod_mcmc.f90 \
        $(SRC_F_DIR)/bayesRCO.f90

# Object files
C_OBJ = $(patsubst $(SRC_C_DIR)/%.c, $(OBJ_DIR)/c/%.o, $(C_SRC))
F_OBJ = $(patsubst $(SRC_F_DIR)/%.f90, $(OBJ_DIR)/fortran/%.o, $(F_SRC))

# Executables
C_EXEC = $(BIN_DIR)/bayesRCO_c
F_EXEC = $(BIN_DIR)/bayesRCO_fortran

# Default target
all: directories $(F_EXEC) $(C_EXEC)

# Create directories
directories:
	@mkdir -p $(OBJ_DIR)/c
	@mkdir -p $(OBJ_DIR)/fortran
	@mkdir -p $(BIN_DIR)

# Fortran build
$(F_EXEC): $(F_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJ_DIR)/fortran/RandomDistributions.o: $(SRC_F_DIR)/RandomDistributions.f90
	$(FC) $(FFLAGS) -J$(OBJ_DIR)/fortran -c -o $@ $<

$(OBJ_DIR)/fortran/mod_types.o: $(SRC_F_DIR)/mod_types.f90
	$(FC) $(FFLAGS) -J$(OBJ_DIR)/fortran -c -o $@ $<

$(OBJ_DIR)/fortran/mod_io.o: $(SRC_F_DIR)/mod_io.f90 $(OBJ_DIR)/fortran/mod_types.o
	$(FC) $(FFLAGS) -J$(OBJ_DIR)/fortran -c -o $@ $<

$(OBJ_DIR)/fortran/mod_mcmc.o: $(SRC_F_DIR)/mod_mcmc.f90 $(OBJ_DIR)/fortran/mod_types.o $(OBJ_DIR)/fortran/RandomDistributions.o
	$(FC) $(FFLAGS) -J$(OBJ_DIR)/fortran -c -o $@ $<

$(OBJ_DIR)/fortran/bayesRCO.o: $(SRC_F_DIR)/bayesRCO.f90 \
                               $(OBJ_DIR)/fortran/mod_types.o \
                               $(OBJ_DIR)/fortran/mod_io.o \
                               $(OBJ_DIR)/fortran/mod_mcmc.o
	$(FC) $(FFLAGS) -J$(OBJ_DIR)/fortran -c -o $@ $<

# C build
$(C_EXEC): $(C_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/c/%.o: $(SRC_C_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean directories
