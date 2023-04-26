CPP=/opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12

GSL_LIB=/opt/homebrew/Cellar/gsl/2.7.1
GSL_DIR=/opt/homebrew/Cellar/gsl/2.7.1

RUNTIME_DIR=./runtime
SAMPLER_DIR=./sampler
BIN_DIR=./bin

PARA_TARGET=gemm-t4-pluss-pro-model-ri-omp
SEQ_TARGET=gemm-t4-pluss-pro-model-ri-omp-seq
TNUM=4
CPPFLAG=-DTHREAD_NUM=$(TNUM) -DCHUNK_SIZE=4 -DDS=8 -DCLS=64 #-DDEBUG

check:
	$(CPP) -std=c++17 -O3 -lpthread -fopenmp -I$(GSL_DIR)/include -I$(RUNTIME_DIR) -c $(SAMPLER_DIR)/$(PARA_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib $(PARA_TARGET).o pluss.o pluss_utils.o -lgsl -lgslcblas -lm -lpthread -fopenmp -o $(BIN_DIR)/$(PARA_TARGET)
	$(CPP) -std=c++17 -O3 -lpthread -fopenmp -I$(GSL_DIR)/include -I$(RUNTIME_DIR) -c $(SAMPLER_DIR)/$(SEQ_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib $(SEQ_TARGET).o pluss.o pluss_utils.o -lgsl -lgslcblas -lm -lpthread -fopenmp -o $(BIN_DIR)/$(SEQ_TARGET)

clean:
	rm $(BIN_DIR)/$(PARA_TARGET)
	rm pluss.o pluss_utils.o $(PARA_TARGET).o
	rm $(BIN_DIR)/$(SEQ_TARGET)
	rm $(SEQ_TARGET).o
	rm ../../output.txt

run: check
#	$(BIN_DIR)/$(PARA_TARGET)
#	$(BIN_DIR)/$(SEQ_TARGET)
	for number in 1 2 3 ; do \
		$(BIN_DIR)/$(PARA_TARGET) >> ../../output.txt ; \
		$(BIN_DIR)/$(SEQ_TARGET) >> ../../output.txt ; \
	done