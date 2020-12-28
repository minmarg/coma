INCLUDES = -I$(top_srcdir)/src/ext -I$(top_srcdir)/src/library $(all_includes)
METASOURCES = AUTO
bin_PROGRAMS = coma
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
coma_SOURCES = coma.cpp
coma_LDADD = $(top_builddir)/src/library/libprobox.a \
	$(top_builddir)/src/ext/libpsl.a
