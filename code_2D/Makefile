# Authors: Gabriel Suau, Remi Pegouret, Lucas Trautmann
#
# Version 0.1.0
# 
# Copyright © 2021 Gabriel Suau
# Copyright © 2021 Remi Pegouret
# Copyright © 2021 Lucas Trautmann
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#########################################################
# 			USAGE				#
#							#
# 	- compilation en mode debug : make debug	#
# 	- compilation en mode optimisé : make release	#
#							#
#########################################################

# Compilateur + flags génériques
CC        = g++
CXX_FLAGS = -std=c++11 -I Eigen/Eigen

# Flags d'optimisation et de debug
OPTIM_FLAGS = -O2 -DNDEBUG
DEBUG_FLAGS = -O0 -g -DDEBUG -pedantic -fbounds-check -fdump-core -pg

# Nom de l'exécutable
PROG = main
# Fichiers sources
SRC = main.cpp DataFile.cpp Mesh.cpp Physics.cpp FiniteVolume.cpp TimeScheme.cpp

.PHONY: release debug clean

# Mode release par défaut
release: CXX_FLAGS += $(OPTIM_FLAGS)
release: $(PROG)

# Mode debug
debug: CXX_FLAGS += $(DEBUG_FLAGS)
debug: $(PROG)

# Compilation + édition de liens
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)


# Supprime l'exécutable, les fichiers binaires (.o), les fichiers
# temporaires de sauvegarde (~), et le fichier de profiling (.out)
clean:
	rm -f *.o *~ gmon.out $(PROG)
