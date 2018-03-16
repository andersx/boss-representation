F2PY=/home/andersx/.local/bin/f2py

# Flags for GCC compilers and MKL
COMPILER_FLAGS = --opt='-O3 -fopenmp -O3 -m64 -march=native' --f90flags='-I${MKLROOT}/include'
LINKER_FLAGS = -lgomp -lpthread -lm -ldl

all: fboss.so

fboss.so: fboss.f90
	${F2PY} -c fboss.f90 -m fboss $(COMPILER_FLAGS) $(LINKER_FLAGS)

clean:
	rm -f fboss.so
