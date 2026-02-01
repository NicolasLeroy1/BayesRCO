# Compiler
FC = gfortran

# Directories
SRC_DIR = new_src
REF_DIR = src
BIN_DIR = bin
BUILD_DIR = build
BUILD_REF_DIR = build_ref

# Executables
TARGET = $(BIN_DIR)/bayesRCO
TARGET_REF = $(BIN_DIR)/bayesRCO_ref

# Source files for the new version
# Order matters slightly for compilation, but dependencies are explicit below
MODULE_FILES = mod_defs.f90 mod_data.f90 mod_random.f90 mod_cmd.f90 mod_io.f90 \
               mod_standardize.f90 mod_stats.f90 mod_mcmc_utils.f90 \
               mod_mcmc_mixture.f90 mod_mcmc_additive.f90 mod_mcmc_bayesCpi.f90 mod_mcmc.f90
MODULES = $(addprefix $(SRC_DIR)/, $(MODULE_FILES))
MAIN = $(SRC_DIR)/bayesRCO.f90
SOURCES = $(MODULES) $(MAIN)
OBJECTS = $(addprefix $(BUILD_DIR)/, $(MODULE_FILES:.f90=.o)) $(BUILD_DIR)/bayesRCO.o

# Source files for the reference version
REF_MODULE_FILES = RandomDistributions.f90 baymodsRCO.f90
REF_MODULES = $(addprefix $(REF_DIR)/, $(REF_MODULE_FILES))
REF_MAIN = $(REF_DIR)/bayesRCO.f90
REF_SOURCES = $(REF_MODULES) $(REF_MAIN)
REF_OBJECTS = $(addprefix $(BUILD_REF_DIR)/, $(REF_MODULE_FILES:.f90=.o)) $(BUILD_REF_DIR)/bayesRCO.o

.PHONY: all clean ref new verify

all: ref new

# New modularized version (Fortran)
new: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(FC) -O3 -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) -O3 -J$(BUILD_DIR) -c -o $@ $<

# Explicit Module Dependencies (from 'use' statements)
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

# Helper to ensure directories exist
$(BIN_DIR) $(BUILD_DIR) $(BUILD_REF_DIR):
	mkdir -p $@

# Verification
verify: all
	./verify_equivalence.sh

clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR) $(BUILD_REF_DIR) out_ref out_new
	find . -name "*.mod" -delete
