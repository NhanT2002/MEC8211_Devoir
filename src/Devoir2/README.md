Il faut avoir la librairie eigen3 d'installer sur la machine et le mettre dans /usr/include/eigen3.
Modifier le fichier parameters.txt avant d'executer.

## Pour compiler :
```bash
g++ -std=c++20 -Wall -Wextra -I/usr/include/eigen3 -fopenmp -O3 -ftree-vectorize -march=native -mtune=native -flto=auto -fexceptions -fno-strict-aliasing -Wno-unused-variable main.cpp -o main
```

## Pour executer le code :
```bash
./main parameters.txt
```
## Pour executer l'étude de convergence :
```bash
chmod +x convergence_analysis.sh
./convergence_analysis.sh
```
