INCLUDES = -I$(top_srcdir)/src/library $(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
noinst_LIBRARIES = libpmpi.a
libpmpi_a_SOURCES = MessageDispatcher.cpp
libpmpi_a_LIBADD = $(top_builddir)/src/library/libprobox.a
