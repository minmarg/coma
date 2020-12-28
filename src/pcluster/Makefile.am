INCLUDES = -I$(top_srcdir)/src/library $(all_includes)
METASOURCES = AUTO
bin_PROGRAMS = pcluster
pcluster_SOURCES = pcluster.cpp
pcluster_LDADD = $(top_builddir)/src/library/libprobox.a \
	$(top_builddir)/src/ext/libpsl.a
