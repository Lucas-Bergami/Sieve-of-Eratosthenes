EXEC = Crivo_de_Eratostenes       # paralelo
SEQ_EXEC = Crivo_de_Eratostenes_seq  # sequencial

SRC = Crivo_de_Eratostenes.c      # paralelo
SEQ_SRC = Crivo_de_Eratostenes_Sequencial.c  # sequencial

MPICC = mpicc
CC = gcc

NP = 4
N = 1000000

# ----------------------
# Versão paralela MPI
# ----------------------
all: $(EXEC)

$(EXEC): $(SRC)
	$(MPICC) -o $(EXEC) $(SRC) -lm
	chmod +x $(EXEC)

run: $(EXEC)
	mpirun -np $(NP) ./$(EXEC) $(N)

profile: $(SRC)
	$(MPICC) -pg -o $(EXEC) $(SRC) -lm
	chmod +x $(EXEC)
	mpirun -np $(NP) ./$(EXEC) $(N)
	@echo "Running gprof on gmon.out. Note: This may only reflect the profile of the last process to finish."
	gprof $(EXEC) gmon.out > gprof_report.txt
	@echo "Relatório de profiling gerado em gprof_report.txt"

# ----------------------
# Versão sequencial
# ----------------------
seq: $(SEQ_EXEC)

$(SEQ_EXEC): $(SEQ_SRC)
	$(CC) -o $(SEQ_EXEC) $(SEQ_SRC) -lm
	chmod +x $(SEQ_EXEC)

seq_run: $(SEQ_EXEC)
	./$(SEQ_EXEC) $(N)

seq_profile: $(SEQ_SRC)
	$(CC) -pg -o $(SEQ_EXEC) $(SEQ_SRC) -lm
	chmod +x $(SEQ_EXEC)
	./$(SEQ_EXEC) $(N)
	gprof $(SEQ_EXEC) gmon.out > gprof_seq_report.txt
	@echo "Relatório de profiling sequencial gerado em gprof_seq_report.txt"

# ----------------------
# Limpeza
# ----------------------
clean:
	rm -f $(EXEC) $(SEQ_EXEC) primos.txt gmon.out gprof_report.txt gprof_seq_report.txt
