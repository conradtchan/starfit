.phoney: install

install:
	rm -rf ./build
	make -C ./src/starfit/fitness clean
	SETUPTOOLS_ENABLE_FEATURES="legacy-editable" pip3 install -e .[testing]
