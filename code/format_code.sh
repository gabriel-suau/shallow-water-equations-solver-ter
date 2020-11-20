#!/bin/bash

# Format the code using Artistic Style.

echo "Formatting your code with clang-format..."
clang-format -style=file -fallback-style=none src/*.h src/*.cpp
