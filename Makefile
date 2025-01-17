# Makefile for building the simulation program

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# Target executable
TARGET = simulation

# Source files
SRC = modver.cpp

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean build artifacts
clean:
	rm -f $(TARGET) *.o results.dat

# Run the simulation
run: $(TARGET)
	./$(TARGET)
