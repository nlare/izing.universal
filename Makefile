all:
	g++ -fopenmp -g main.cpp -o main -ffast-math -flto -march=native -msse4.2 -xhost -O3
clean:
	rm main APPR.DAT APPR.SORT.DAT RELAX.DAT
