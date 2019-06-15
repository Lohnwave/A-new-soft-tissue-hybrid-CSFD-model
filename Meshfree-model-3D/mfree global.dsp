# Microsoft Developer Studio Project File - Name="mfree global" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=mfree global - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mfree global.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mfree global.mak" CFG="mfree global - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mfree global - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "mfree global - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "mfree global - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "mfree global - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "mfree global - Win32 Release"
# Name "mfree global - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\bandsolver.f90
# End Source File
# Begin Source File

SOURCE=.\cellgausspoints.f90
DEP_F90_CELLG=\
	".\parameter.h"\
	
# End Source File
# Begin Source File

SOURCE=".\cmp-radial-basis.f90"
# End Source File
# Begin Source File

SOURCE=.\essentialbc.f90
DEP_F90_ESSEN=\
	".\parameter.h"\
	
# End Source File
# Begin Source File

SOURCE=.\gausscoefficient.f90
# End Source File
# Begin Source File

SOURCE=.\gausssolver.f90
# End Source File
# Begin Source File

SOURCE=.\getdisplacement.f90
# End Source File
# Begin Source File

SOURCE=.\getinvasy.f90
# End Source File
# Begin Source File

SOURCE=.\getstress.f90
# End Source File
# Begin Source File

SOURCE=.\input.f90
# End Source File
# Begin Source File

SOURCE=".\mfree global.f90"
DEP_F90_MFREE=\
	".\parameter.h"\
	".\variables.h"\
	
# End Source File
# Begin Source File

SOURCE=.\naturalbcconcentrated.f90
# End Source File
# Begin Source File

SOURCE=.\naturalbcdistributed.f90
DEP_F90_NATUR=\
	".\parameter.h"\
	
# End Source File
# Begin Source File

SOURCE=.\pointstiffnessmatrix.f90
# End Source File
# Begin Source File

SOURCE=.\sem.f90
NODEP_F90_SEM_F=\
	".\Debug\m_gauss.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\seme.f90
# End Source File
# Begin Source File

SOURCE=.\shapefunction.f90
# End Source File
# Begin Source File

SOURCE=.\solverband.f90
# End Source File
# Begin Source File

SOURCE=.\supportdomain.f90
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\parameter.h
# End Source File
# Begin Source File

SOURCE=.\variables.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\Input175.dat
# End Source File
# Begin Source File

SOURCE=.\result.dat
# End Source File
# End Group
# End Target
# End Project
