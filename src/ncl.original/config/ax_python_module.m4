dnl @synopsis AC_PYTHON_MODULE(modname[, fatal[, msg]])
dnl
dnl Checks for Python module.
dnl
dnl If fatal is non-empty then absence of a module will trigger an
dnl error, if fatal is empty then PYTHON will be set to an empty string
dnl
dnl If errormsg is non-empty, the message will be printed
dnl instead of a generic error message.
dnl
dnl @category InstalledPackages
dnl @author Andrew Collier <colliera@nu.ac.za>.
dnl @version 2004-07-14
dnl @license AllPermissive

AC_DEFUN([AC_PYTHON_MODULE],[
	if test -z "PYTHON"
	then
		AC_PATH_PROG([PYTHON], [python], [], [$PATH])
	fi
	AC_MSG_CHECKING(python module: $1)
	if test -z "$PYTHON"
	then
		false
	else
		"$PYTHON" -c "import $1" 2>/dev/null
	fi
	if test $? -eq 0;
	then
		AC_MSG_RESULT(yes)
		eval AS_TR_CPP(HAVE_PYMOD_$1)=yes
	else
		AC_MSG_RESULT(no)
		eval AS_TR_CPP(HAVE_PYMOD_$1)=no
		
		if test -n "$2"
		then
			if test -n "$3"
			then
				AC_MSG_ERROR($3)
			else
				AC_MSG_ERROR(failed to find required module $1)
			fi
			exit 1
		else
			PYTHON=""
		fi
	fi
])
