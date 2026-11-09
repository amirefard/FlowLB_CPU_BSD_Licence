# Compiler
CXX = mpicxx

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Target
TARGET = $(BIN_DIR)/FlowLB_CPU

# Flags
CXXFLAGS = -O2 -Wall -std=c++17 -Iinclude -MMD -MP

# Sources
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# Objects
OBJ = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compile
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Dependencies
-include $(OBJ:.o=.d)