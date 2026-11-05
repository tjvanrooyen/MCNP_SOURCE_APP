# Make file

# Do not edit.

dummy : \
    ../exe/mcnpsource
	@echo Targets updated.

# New target

../exe/mcnpsource : \
    obj/m_KINDS.o \
    obj/FileHandleMod.o \
    obj/mcnpsource.o
	mkexe ../exe/mcnpsource \
    obj/m_KINDS.o \
    obj/FileHandleMod.o \
    obj/mcnpsource.o

obj/m_KINDS.o : m_KINDS.f90
	cmplf90 m_KINDS

obj/FileHandleMod.o : FileHandleMod.f90
	cmplf90 FileHandleMod

obj/mcnpsource.o : mcnpsource.f90
	cmplf90 mcnpsource

m_KINDS.f90 \
FileHandleMod.f90 \
mcnpsource.f90 : mcnpsource.src
	clextsrc mcnpsource.src
