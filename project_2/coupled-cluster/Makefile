RM=rm -rf

BDIR=build
COULOMBDIR=coupled_cluster/matrix_elements
CCDIR=coupled_cluster/schemes

build:
	python3 setup.py build_ext

install:
	@pip3 install . --upgrade

installr:
	@pip3 install . --upgrade -r requirements.txt

.PHONY: clean

clean:
	$(RM) $(BDIR) $(COULOMBDIR)/*.cpp $(CCDIR)/*.cpp .cache
