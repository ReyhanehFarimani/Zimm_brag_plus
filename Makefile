CXX      ?= g++
CXXFLAGS ?= -std=c++11 -O3 -march=native -Wall -Wextra
SRC      := $(wildcard src/*.cpp)
OBJ      := $(SRC:src/%.cpp=build/%.o)
TARGET   := zimm

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

build:
	mkdir -p build

run: $(TARGET)
	./$(TARGET) input.dat

clean:
	rm -rf build $(TARGET)

-include $(OBJ:.o=.d)
.PHONY: all run clean
