F2PY=/home/andersx/.local/bin/f2py

all: fboss.so

fboss.so: fboss.f90
	${F2PY} -c fboss.f90 -m fboss
