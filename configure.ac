#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([tulip], [2], [krister.henriksson@helsinki.fi])
AC_CONFIG_SRCDIR([src/param-pot.cpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CXX
# AC_PROG_CC


# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([rt], [clock_gettime])
AC_CHECK_LIB([symspg], [spg_get_dataset])
AC_CHECK_LIB([utils], [aborterror])
AC_CHECK_LIB([omp], [omp_get_num_threads])

# Checks for header files.

AC_C_INLINE
AC_TYPE_SIZE_T


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
