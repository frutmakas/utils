# Microsoft Developer Studio Project File - Name="utils" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=utils - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "utils.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "utils.mak" CFG="utils - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "utils - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "utils - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE "utils - Win32 200x400" (based on "Win32 (x86) Console Application")
!MESSAGE "utils - Win32 Release_2" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "utils - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /G6 /MT /W3 /GX /O2 /Op /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT CPP /Fr
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# SUBTRACT LINK32 /profile /incremental:yes

!ELSEIF  "$(CFG)" == "utils - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /profile /debug /machine:I386 /out:"Debug/stfbc.exe"

!ELSEIF  "$(CFG)" == "utils - Win32 200x400"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "utils___Win32_200x400"
# PROP BASE Intermediate_Dir "utils___Win32_200x400"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "utils___Win32_200x400"
# PROP Intermediate_Dir "utils___Win32_200x400"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"200x400/utils.exe"
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"200x400/utils.exe"

!ELSEIF  "$(CFG)" == "utils - Win32 Release_2"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "utils___Win32_Release_2"
# PROP BASE Intermediate_Dir "utils___Win32_Release_2"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "utils___Win32_Release_2"
# PROP Intermediate_Dir "utils___Win32_Release_2"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /G6 /MT /W3 /GX /O2 /Op /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /G6 /MT /W3 /GX /O2 /Op /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT CPP /Fr
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# SUBTRACT BASE LINK32 /profile /incremental:yes
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# SUBTRACT LINK32 /profile /incremental:yes

!ENDIF 

# Begin Target

# Name "utils - Win32 Release"
# Name "utils - Win32 Debug"
# Name "utils - Win32 200x400"
# Name "utils - Win32 Release_2"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\stbc\conversion.cpp
# End Source File
# Begin Source File

SOURCE=..\tools\ext_utilitis.cpp
# End Source File
# Begin Source File

SOURCE=..\fading\fadingtdma.cpp
# End Source File
# Begin Source File

SOURCE=..\fft\fft.cpp
# End Source File
# Begin Source File

SOURCE=..\queue\fifo.cpp
# End Source File
# Begin Source File

SOURCE=..\cdma\gold.cpp
# End Source File
# Begin Source File

SOURCE=..\sort\heapsort.cpp
# End Source File
# Begin Source File

SOURCE=..\tools\interleaver.cpp
# End Source File
# Begin Source File

SOURCE=..\interpol\interpol.cpp
# End Source File
# Begin Source File

SOURCE=..\coding\block\ldpc.cpp
# End Source File
# Begin Source File

SOURCE=..\queue\list.cpp
# End Source File
# Begin Source File

SOURCE=..\nr\nrutil.cpp
# End Source File
# Begin Source File

SOURCE=..\modulation\ofdm.cpp
# End Source File
# Begin Source File

SOURCE=..\tools\pilot.cpp
# End Source File
# Begin Source File

SOURCE=..\modulation\psk.cpp
# End Source File
# Begin Source File

SOURCE=..\sort\quicksort.cpp
# End Source File
# Begin Source File

SOURCE=..\rand\randgen.cpp
# End Source File
# Begin Source File

SOURCE=..\cdma\spread.cpp
# End Source File
# Begin Source File

SOURCE=..\tools\tools.cpp
# End Source File
# Begin Source File

SOURCE=..\tools\utilitis.cpp
# End Source File
# Begin Source File

SOURCE=..\cdma\walsh.cpp
# End Source File
# Begin Source File

SOURCE=..\coding\zigzag.cpp
# End Source File
# Begin Source File

SOURCE=..\coding\zigzag.h
# End Source File
# Begin Source File

SOURCE=.\zigzag_test.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\tools\all.h
# End Source File
# Begin Source File

SOURCE=.\all.h
# End Source File
# Begin Source File

SOURCE=..\tree\btree.h
# End Source File
# Begin Source File

SOURCE=..\stbc\conversion.h
# End Source File
# Begin Source File

SOURCE=..\sort\exception.h
# End Source File
# Begin Source File

SOURCE=..\tools\ext_utilitis.h
# End Source File
# Begin Source File

SOURCE=..\fading\fadingtdma.h
# End Source File
# Begin Source File

SOURCE=..\fft\fft.h
# End Source File
# Begin Source File

SOURCE=..\queue\fifo.h
# End Source File
# Begin Source File

SOURCE=..\globaldef.h
# End Source File
# Begin Source File

SOURCE=..\cdma\gold.h
# End Source File
# Begin Source File

SOURCE=..\sort\heapsort.h
# End Source File
# Begin Source File

SOURCE=..\tools\interleaver.h
# End Source File
# Begin Source File

SOURCE=..\interpol\interpol.h
# End Source File
# Begin Source File

SOURCE=..\coding\block\ldpc.h
# End Source File
# Begin Source File

SOURCE=..\queue\list.h
# End Source File
# Begin Source File

SOURCE=..\nr\nrutil.h
# End Source File
# Begin Source File

SOURCE=..\modulation\ofdm.h
# End Source File
# Begin Source File

SOURCE=..\tools\pilot.h
# End Source File
# Begin Source File

SOURCE=..\modulation\psk.h
# End Source File
# Begin Source File

SOURCE=..\sort\quicksort.h
# End Source File
# Begin Source File

SOURCE=..\rand\randgen.h
# End Source File
# Begin Source File

SOURCE=..\cdma\spread.h
# End Source File
# Begin Source File

SOURCE=..\tools\superclass.h
# End Source File
# Begin Source File

SOURCE=..\tools\tools.h
# End Source File
# Begin Source File

SOURCE=..\tools\utilitis.h
# End Source File
# Begin Source File

SOURCE=..\cdma\walsh.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
