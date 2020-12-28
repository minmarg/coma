INCLUDES = -I$(top_srcdir)/src/library $(all_includes)
METASOURCES = AUTO
bin_PROGRAMS = makepro
makepro_SOURCES = makepro.cpp
makepro_LDADD = $(top_builddir)/src/library/libprobox.a \
	$(top_builddir)/src/ext/libpsl.a
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
