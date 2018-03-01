RM=rm -rf

BDIR=build
CYTHONDIR=vmc

build:
	@python3 setup.py build_ext

install:
	@pip install . --upgrade

.PHONY: clean

clean:
	$(RM) $(BDIR) $(CYTHONDIR)/*.cpp .cache
