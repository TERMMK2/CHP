# Compilateur utilisé
CC=g++
CC_para=mpic++ 


run : mainPara.cc fonction.cpp Laplacian2DPara.cpp fonction.h Laplacian2DPara.h
	$(CC_para) -std=c++11  mainPara.cc fonction.cpp Laplacian2DPara.cpp  -o run

#si on a des trucs a tester :
test : test.cc 
	$(CC) test.cc  $(CXX_FLAGS) -o run_test




# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ run
