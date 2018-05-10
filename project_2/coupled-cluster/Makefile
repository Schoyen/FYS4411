RM=rm -rf

BDIR=build
CYTHONDIR=coupled_cluster/matrix_elements

build:
	python3 setup.py build_ext

install:
	@pip3 install . --upgrade

installr:
	@pip3 install . --upgrade -r requirements.txt

.PHONY: clean

clean:
	$(RM) $(BDIR) $(CYTHONDIR)/*.cpp .cache
