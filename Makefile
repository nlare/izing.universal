all:
	g++ -fopenmp -g main.cpp -o main
clean:
	rm main APPR.DAT APPR.SORT.DAT RELAX.DAT