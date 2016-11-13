SUBDIRS = canal fading sort queue cdma modulation rand coding fft interpol nr stbc map nspmatrice tools utils

all:
	if test -e test_utilitis; then rm test_utilitis;fi
	if test -e test_ldpc; then rm test_ldpc;fi
	if test -e emtest; then rm emtest;fi
	if test -e em2; then rm em2;fi
	if test -e mimostbc; then rm mimostbc; fi
	if test -e mccdma2; then rm mccdma2; fi
	if test -e zigzag_test.exe; then rm zigzag_test.exe; fi
	list='$(SUBDIRS)'; for subdir in $$list; do \
	echo "Making in $$subdir"; \
	(cd $$subdir && if test -e Makefile; then \
	make; \
	else echo "Nothing to do in $$subdir"; fi)\
	done; \

debug: 
	if test -e test_utilitis; then rm test_utilitis;fi
	if test -e test_ldpc; then rm test_ldpc;fi
	if test -e emtest; then rm emtest;fi
	if test -e em2; then rm em2;fi
	list='$(SUBDIRS)'; for subdir in $$list; do \
	echo "Making in $$subdir"; \
	(cd $$subdir && if test -e Makefile.debug; then \
	make -f Makefile.debug; \
	else echo "Nothing to do in $$subdir"; fi)\
	done; \

clean:
	for i in `find . -name "*.o"`; do (rm -fv $$i); done
	if test -e test_utilitis; then rm test_utilitis;fi
	if test -e test_ldpc; then rm test_ldpc;fi
	if test -e emtest; then rm emtest;fi
	if test -e em2; then rm em2;fi
	if test -e mimostbc; then rm mimostbc; fi
	if test -e zigzag_test.exe; then rm zigzag_test.exe; fi
