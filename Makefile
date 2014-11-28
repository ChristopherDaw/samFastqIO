all:
	$(MAKE) -C src
	mkdir -p bin
	mv src/sfio bin/sfio

debug:
	$(MAKE) -C src debug
	mkdir -p bin
	mv src/sfio bin/sfio

clean:
	$(MAKE) -C src clean
	rm -f bin/sfio
