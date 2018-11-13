DOXYGEN:=$(shell which doxygen)
GFORTRAN:=$(shell which gfortran)

all: shan fortran

fortran:
	$(MAKE) -C src fortran

shan:
	$(MAKE) -C src shan

docs:
	@if test "$(DOXYGEN)" = ""; then \
		echo "Doxygen not found."; \
		echo "Install doxygen to be able to generate documentation."; \
		false; \
	fi
	doxygen Doxyfile

clean:
	$(MAKE) -C src clean

.PHONY: all tests docs clean 
