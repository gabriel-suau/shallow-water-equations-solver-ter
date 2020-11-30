#!/bin/bash

# Format the code using clang-format.

echo "Formatting your code with clang-format..."
clang-format -style=file -fallback-style=none *.h *.cpp
