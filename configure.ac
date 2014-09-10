#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)

AC_INIT([src/sax/sax.h])
AM_INIT_AUTOMAKE([squad],[2.0],[chandola@cs.umn.edu])
AM_CONFIG_HEADER(config.h)
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_CHECK_TOOL(RANLIB, ranlib)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_VPRINTF

AC_CHECK_FUNC(round, ,AC_CHECK_LIB(m,round))
AC_OUTPUT([Makefile
	src/Makefile
        src/sax/Makefile
        src/ar/Makefile
        src/box/Makefile
        src/clusterc/Makefile
        src/clusterd/Makefile
        src/data-generator/Makefile
        src/fsa/Makefile
        src/hmm/Makefile
        src/knnc/Makefile
        src/knnd/Makefile
        src/svr/Makefile
        src/ocsvm/Makefile
        src/winc/Makefile
        src/wind/Makefile]
)
