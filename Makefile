.phoney: install

install:
	export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
	rm -rf ./build
	make -C ./src/starfit/fitness clean
	pip3 install -e .[testing]
