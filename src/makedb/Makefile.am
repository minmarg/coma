INCLUDES = -I$(top_srcdir)/src/library $(all_includes)
METASOURCES = AUTO
bin_PROGRAMS = makedb
makedb_SOURCES = makedb.cpp
makedb_LDADD = $(top_builddir)/src/library/libprobox.a \
	$(top_builddir)/src/ext/libpsl.a
