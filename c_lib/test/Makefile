CPP=g++

GSL_LIB=/usr/include/gsl
GSL_DIR=/usr/include/gsl

RUNTIME_DIR=./runtime
SAMPLER_DIR=./sampler
BIN_DIR=./bin

PARA_TARGET=gemm-t4-pluss-pro-model-ri-omp
SEQ_TARGET=gemm-t4-pluss-pro-model-ri-omp-seq
OPT_TARGET=gemm-t4-pluss-pro-model-ri-opt
R10_TARGET=gemm-t4-pluss-pro-model-rs-ri-opt-r10
TNUM=4
CPPFLAG=-DTHREAD_NUM=$(TNUM) -DCHUNK_SIZE=4 -DDS=8 -DCLS=64 #-DDEBUG

check:
	$(CPP) -std=c++17 -O3 -flto=auto -fvisibility=hidden -lpthread -fopenmp -I$(GSL_DIR)/include -I$(RUNTIME_DIR) -c $(SAMPLER_DIR)/$(PARA_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib -flto=auto -fvisibility=hidden $(PARA_TARGET).o pluss.o pluss_utils.o -lgsl -lgslcblas -lm -lpthread -fopenmp -o $(BIN_DIR)/$(PARA_TARGET)
	$(CPP) -std=c++17 -O3 -flto=auto -fvisibility=hidden -lpthread -fopenmp -I$(GSL_DIR)/include -I$(RUNTIME_DIR) -c $(SAMPLER_DIR)/$(SEQ_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib -flto=auto -fvisibility=hidden $(SEQ_TARGET).o pluss.o pluss_utils.o -lgsl -lgslcblas -lm -lpthread -fopenmp -o $(BIN_DIR)/$(SEQ_TARGET)
	$(CPP) -std=c++17 -flto=auto -fvisibility=hidden -O3 -g3 -fno-omit-frame-pointer -lpthread -fopenmp -march=native -I/localdisk/local/include -I$(RUNTIME_DIR) -I../../parallel-hashmap -c $(SAMPLER_DIR)/$(OPT_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils_v2.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib -flto=auto -fvisibility=hidden $(OPT_TARGET).o pluss.o pluss_utils_v2.o -lgsl -lgslcblas -lm -lpthread -fopenmp -o $(BIN_DIR)/$(OPT_TARGET)
	$(CPP) -std=c++17 -flto=auto -fvisibility=hidden -fopenmp -lpthread -O3 -march=native -I$(GSL_DIR)/include -I$(RUNTIME_DIR) -c $(SAMPLER_DIR)/$(R10_TARGET).cpp $(RUNTIME_DIR)/pluss.cpp $(RUNTIME_DIR)/pluss_utils.cpp $(CPPFLAG)
	$(CPP) -L$(GSL_LIB)/lib -flto=auto -fvisibility=hidden -fopenmp -lpthread $(R10_TARGET).o pluss.o pluss_utils.o -lgsl -lgslcblas -lm -o $(BIN_DIR)/$(R10_TARGET)

clean:
	rm $(BIN_DIR)/$(PARA_TARGET)
	rm *.o
	rm $(BIN_DIR)/$(SEQ_TARGET) $(BIN_DIR)/$(OPT_TARGET)
	rm ../../output.txt
	rm runtime/*.o

speed: check
	$(BIN_DIR)/$(SEQ_TARGET) >> ../../output.txt speed
	$(BIN_DIR)/$(PARA_TARGET) >> ../../output.txt speed
	$(BIN_DIR)/$(OPT_TARGET) >> ../../output.txt speed

acc: check
	$(BIN_DIR)/$(PARA_TARGET) >> ../../output.txt acc
	$(BIN_DIR)/$(SEQ_TARGET) >> ../../output.txt acc

sample: check
	$(BIN_DIR)/$(R10_TARGET)