CXX ?= g++
CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic

BUILD_DIR := build
TARGET := $(BUILD_DIR)/tsp_solver
SOURCES := src/mst.cpp src/simplex.cpp src/tsp_solver.cpp
HEADERS := include/mst.h include/simplex.h

.PHONY: all clean run smoke-test

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(SOURCES) $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

run: $(TARGET)
	$(TARGET)

smoke-test: $(TARGET)
	printf 'examples/usa-demo/adjacency_matrix.txt\nexamples/usa-demo\ngraph.tex\nusa.tex\nEXIT\n' | $(TARGET)

clean:
	rm -rf $(BUILD_DIR) lp.txt
