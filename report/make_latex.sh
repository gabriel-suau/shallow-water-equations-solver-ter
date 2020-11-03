#!/bin/bash

echo "LaTeX compilation n°1..."
xelatex -shell-escape main.tex
echo "Building glossary..."
makeglossaries main
echo "LaTeX compilation n°2..."
xelatex -shell-escape main.tex
echo "Building bibliography..."
biber main
echo "Building nomenclature..."
makeindex main 
echo "LaTeX compilation n°3..."
xelatex -shell-escape main.tex
echo "LaTeX compilation n°4..."
xelatex -shell-escape main.tex
