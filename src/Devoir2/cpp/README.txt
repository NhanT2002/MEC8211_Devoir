Il faut avoir la librairie eigen3 d'installer sur la machine

Pour compiler :

g++ -std=c++20 -Wall -Wextra -fopenmp -O3 -ftree-vectorize -march=native -mtune=native -flto=auto -fexceptions -fno-strict-aliasing -I/usr/include/eigen3 -Wno-unused-variable main.cpp -o main

Pour executer l'étude de convergence :
chmod +x convergence_analysis.sh
./convergence_analysis.sh
