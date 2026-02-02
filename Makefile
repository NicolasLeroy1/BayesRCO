# Compilers
FC = gfortran
CC = gcc

# Directories
SRC_DIR = new_src
REF_DIR = src
C_SRC_DIR = new_src_c
BIN_DIR = bin
BUILD_DIR = build
BUILD_REF_DIR = $(BUILD_DIR)/ref
BUILD_C_DIR = $(BUILD_DIR)/c
OUT_DIR = outputs
PLOT_DIR = plots

# Include modular definitions
include $(SRC_DIR)/makefile.mk
include $(REF_DIR)/makefile.mk
# include $(C_SRC_DIR)/makefile.mk # No longer exists

# Executables
TARGET = $(BIN_DIR)/bayesRCO
TARGET_REF = $(BIN_DIR)/bayesRCO_ref
TARGET_C = $(BIN_DIR)/bayesRCO_c

# Objects (Derived from modular includes)
OBJECTS = $(addprefix $(BUILD_DIR)/, $(NEW_MODULE_FILES:.f90=.o)) $(BUILD_DIR)/$(NEW_MAIN:.f90=.o)
REF_OBJECTS = $(addprefix $(BUILD_REF_DIR)/, $(REF_MODULE_FILES:.f90=.o)) $(BUILD_REF_DIR)/$(REF_MAIN:.f90=.o)
# C_OBJECTS = $(addprefix $(BUILD_C_DIR)/, $(C_SOURCE_FILES:.c=.o)) # Managed by sub-make

.PHONY: all clean ref new c verify directories

all: directories ref new c

directories: $(BIN_DIR) $(BUILD_DIR) $(BUILD_REF_DIR) $(BUILD_C_DIR) $(OUT_DIR) $(PLOT_DIR)

# New modularized version (Fortran)
new: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(FC) -O3 -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) -O3 -J$(BUILD_DIR) -c -o $@ $<

# Explicit Module Dependencies (Fortran)
$(BUILD_DIR)/mod_data.o: $(BUILD_DIR)/mod_defs.o
$(BUILD_DIR)/mod_random.o: $(BUILD_DIR)/mod_defs.o
$(BUILD_DIR)/mod_cmd.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o
$(BUILD_DIR)/mod_io.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o
$(BUILD_DIR)/mod_standardize.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o
$(BUILD_DIR)/mod_stats.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o
$(BUILD_DIR)/mod_mcmc_utils.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o $(BUILD_DIR)/mod_io.o
$(BUILD_DIR)/mod_mcmc_mixture.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o $(BUILD_DIR)/mod_random.o $(BUILD_DIR)/mod_io.o $(BUILD_DIR)/mod_stats.o $(BUILD_DIR)/mod_mcmc_utils.o
$(BUILD_DIR)/mod_mcmc_additive.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o $(BUILD_DIR)/mod_random.o $(BUILD_DIR)/mod_io.o $(BUILD_DIR)/mod_stats.o $(BUILD_DIR)/mod_mcmc_utils.o
$(BUILD_DIR)/mod_mcmc_bayesCpi.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o $(BUILD_DIR)/mod_random.o $(BUILD_DIR)/mod_io.o $(BUILD_DIR)/mod_stats.o $(BUILD_DIR)/mod_mcmc_utils.o
$(BUILD_DIR)/mod_mcmc.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_mcmc_mixture.o $(BUILD_DIR)/mod_mcmc_additive.o $(BUILD_DIR)/mod_mcmc_bayesCpi.o
$(BUILD_DIR)/bayesRCO.o: $(BUILD_DIR)/mod_defs.o $(BUILD_DIR)/mod_data.o $(BUILD_DIR)/mod_cmd.o $(BUILD_DIR)/mod_io.o $(BUILD_DIR)/mod_stats.o $(BUILD_DIR)/mod_mcmc.o $(BUILD_DIR)/mod_standardize.o

# Reference version
ref: $(TARGET_REF)

$(TARGET_REF): $(REF_OBJECTS) | $(BIN_DIR)
	$(FC) -O3 -o $@ $^

$(BUILD_REF_DIR)/%.o: $(REF_DIR)/%.f90 | $(BUILD_REF_DIR)
	$(FC) -O3 -J$(BUILD_REF_DIR) -c -o $@ $<

# Reference Explicit Dependencies
$(BUILD_REF_DIR)/baymodsRCO.o: $(REF_DIR)/baymodsRCO.f90 $(BUILD_REF_DIR)/RandomDistributions.o
$(BUILD_REF_DIR)/bayesRCO.o: $(REF_DIR)/bayesRCO.f90 $(BUILD_REF_DIR)/RandomDistributions.o $(BUILD_REF_DIR)/baymodsRCO.o

# C version
c: $(TARGET_C)

$(TARGET_C): | $(BIN_DIR)
	$(MAKE) -C $(C_SRC_DIR)
	cp $(C_SRC_DIR)/clibayesrco/clibayesrco $@

# $(BUILD_C_DIR)/%.o: $(C_SRC_DIR)/%.c | $(BUILD_C_DIR)
# 	$(CC) -O3 -I$(C_SRC_DIR) -c -o $@ $<

# Helper to ensure directories exist
$(BIN_DIR) $(BUILD_DIR) $(BUILD_REF_DIR) $(BUILD_C_DIR) $(OUT_DIR) $(PLOT_DIR):
	mkdir -p $@

# Verification
verify: all
	./tests/equivalence/run_test.sh

clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR) $(OUT_DIR) $(PLOT_DIR) out_ref out_new out_c out_benchmark debug_run debug_c
	find . -name "*.mod" -delete
	find . -name "*.o" -delete

