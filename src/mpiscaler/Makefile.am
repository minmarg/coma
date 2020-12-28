INCLUDES = -I$(top_srcdir)/src/library -I$(top_srcdir)/src/lmpi $(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
bin_SCRIPTS = mpiscaler.sh
bin_PROGRAMS = mpiscaler
mpiscaler_SOURCES = PScaler.cpp ParallelUniversalScoreMatrix.cpp mpiscaler.cpp
mpiscaler_LDADD = $(top_builddir)/src/library/libprobox.a \
	$(top_builddir)/src/lmpi/libpmpi.a $(top_builddir)/src/ext/libpsl.a -lmpich -lpthread
