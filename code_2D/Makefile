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
