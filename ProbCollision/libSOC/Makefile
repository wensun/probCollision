PACKAGE = libSOC.10.1

all: swift lib tests doc

lib::
	make -C lib

tests::
	-make -C test/array
	-make -C test/ors_editor
	-make -C test/ors_fileConverter
	-make -C test/soc
	-make -C test/09-rlss
	-make -C test/10-dai

swift::
	make -C extern/SWIFT++_1.2/src

doc::
	make -C doc guide doxy

clean::
	make -C lib clean
	make -C extern/SWIFT++_1.2/src clean
	make -C test/array clean
	make -C test/ors_editor clean
	make -C test/ors_fileConverter clean
	make -C test/soc clean
	make -C test/09-rlss clean
	make -C test/10-dai clean
	make -C doc clean

zip::
	cd ..;  rm -f $(PACKAGE).tgz;  tar cvzf $(PACKAGE).tgz libSOC --exclude "*.svn"
