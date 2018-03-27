# Compilateur utilisé
CC=g++
CC_para=mpic++

#Chemin Eigen
path=~/Tp\ info/2A/EigenLibrary/Eigen

# Options en mode optimisé - La variable NDEBUG est définie
OPTIM_FLAG = -O3 -DNDEBUG -I $(path) -std=c++11 -Wall
# Options en mode debug
DEBUG_FLAG = -g -I $(path) -std=c++11 -Wall

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.cc Laplacian2D.cpp Datafile.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)

#si on a des trucs a tester :
test : test.cc 
	$(CC) test.cc  $(CXX_FLAGS) -o run_test

CG : main.cc
	$(CC_para) main.cc fonction.cpp -o run_CG


# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
