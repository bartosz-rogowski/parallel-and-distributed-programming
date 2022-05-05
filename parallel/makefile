CC=mpicc
CFLAGS=
DEPS=
OBJS=main.o
OUT=main
SHELL := /bin/bash


$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)
	
.c.o:
	$(CC) -c -o $@ $< $(CFLAGS) 

prepare:
	@echo  "przygotowywanie środowiska..."
	@source /opt/nfs/config/source_mpich32.sh
	@/opt/nfs/config/station206_name_list.sh 1 16 > nodes
#	 @echo -e "\033[1;32mliczba dostępnych węzłów: $(grep "" -c nodes)\033[0m"
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	rm -f *.txt
	

mem:
	valgrind ./$(OUT)

run:
	mpiexec -f nodes -n $(n) ./$(OUT) $(in) $(out)
run01:
	mpiexec -f nodes -n 5 ./$(OUT) example/data01.txt mst_res.txt

run02:
	mpiexec -f nodes -n 10 ./$(OUT) example/data02.txt mst_res.txt

run03:
	mpiexec -f nodes -n 10 ./$(OUT) example/data03.txt mst_res.txt

run04:
	mpiexec -f nodes -n 10 ./$(OUT) example/data04.txt mst_res.txt

help:
	@echo -e "\033[1;33m============= makefile help =============\033[0m"
	@echo -e "\033[1;32mmake\033[0m       - kompiluje i linkuje pliki"
	@echo -e "\033[1;32mmake run\033[0m   - uruchamia program (plik wykonywalny)"
	@echo -e "\033[1;32mmake mem\033[0m   - sprawdza wycieki pamięci za pomocą valgrinda"
	@echo -e "\033[1;32mmake clean\033[0m - usuwa pliki .o oraz plik wykonywalny"
	@echo -e "\033[1;32mmake help\033[0m  - wyświetla pomoc"
	@echo -e "\033[1;33m=========================================\033[0m"