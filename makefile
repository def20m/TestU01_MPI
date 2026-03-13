CC=g++
DIR=${HOME}

seqcrush: exp extest.cpp
	$(CC) -std=c++14 -Ofast extest.cpp -o SeqCrush -I$(DIR)/include -I$(DIR)/sprng5/include -L$(DIR)/lib  -L$(DIR)/sprng5/lib -lmylib -ltestu01 -lsprng
mpicrush: exp ompitest.cpp
	mpic++ -Ofast ompitest.cpp -o MPIcrush -I$(DIR)/include -I$(DIR)/sprng5/include -L$(DIR)/lib  -L$(DIR)/sprng5/lib -ltestu01 -lsprng
test67: exp test67.cpp
	$(CC) -std=c++14 -Ofast test67.cpp -o Test67 -I$(DIR)/include -I$(DIR)/sprng5/include -L$(DIR)/lib  -L$(DIR)/sprng5/lib -lmylib -ltestu01 -lsprng
clean:
	rm -f SeqCrush
	rm -f MPIcrush
	rm -f Test67
exp:
	bash export.sh
