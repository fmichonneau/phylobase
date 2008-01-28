include ${R_HOME}/etc${R_ARCH}/Makeconf
SRC = Rcpp.cpp
OBJ = Rcpp.o
LIB = libRcpp.a

.PHONY: all

all: $(LIB)

clean:
	rm -f $(OBJ) $(LIB)

$(OBJ): $(SRC)

$(LIB): $(OBJ)
	ar crs $(LIB) $(OBJ)
