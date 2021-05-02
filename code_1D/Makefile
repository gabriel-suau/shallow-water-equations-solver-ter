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

# Verbosity level (0,1,2)
# 	0 = Beginning, error and ending logs (not verbose)
# 	1 = All of the above + saving logs (more verbose)
VERBOSITY_LEVEL = 0

CXX_FLAGS += -DVERBOSITY=$(VERBOSITY_LEVEL)

# Flags d'optimisation et de debug
OPTIM_FLAGS = -O2 -DNDEBUG
DEBUG_FLAGS = -O0 -g -DDEBUG -pedantic

# Nom de l'exécutable
PROG = main
# Fichiers sources
SRC = main.cpp DataFile.cpp Mesh.cpp Physics.cpp FiniteVolume.cpp TimeScheme.cpp

# Mode release par défaut
.PHONY: release
release: CXX_FLAGS += $(OPTIM_FLAGS)
release: $(PROG)

# Mode debug
.PHONY: debug
debug: CXX_FLAGS += $(DEBUG_FLAGS)
debug: $(PROG)

# Compilation + édition de liens
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
