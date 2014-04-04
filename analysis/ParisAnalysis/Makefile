
# root definitions 
#
include Makefile.arch
#
# modifications and add ons
#
SrcSuf = C
CXXFLAGS += -Ianalysis

#
#
THEEVENTSLIB    = libParisEvents.$(DllSuf)
THEANALYSISLIB  = libParisAnalysis.$(DllSuf)

ROOTEVENTSO :=  ParisEventsDict.o ParisEvents.o 
ROOTANALYSISSO := ID.o Algo.o ParisAlgo.o ParisEventReader.o Analyser.o

#ID.o Algo.o ParisAlgo.o
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

#
# main target
#
#  the library
#
all: $(THEEVENTSLIB)  
#all: $(THEEVENTSLIB) $(THEANALYSISLIB) 

#	
# rules to build this library
#
#  target to build the library ... maybe better to use MAKELIB ..
#
ifeq ($(ARCH),win32gcc)
EXPLLINKLIBS += 
endif
$(THEEVENTSLIB): $(ROOTEVENTSO)
	@echo " ---> Building $@ ... " 
ifeq ($(ARCH),aix)
		@/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		@/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
          $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif 
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "  ---> $@ done"

#	
# rules to build this library
#
#  target to build the library ... maybe better to use MAKELIB ..
#
ifeq ($(ARCH),win32gcc)
EXPLLINKLIBS += 
endif
$(THEANALYSISLIB): $(ROOTANALYSISSO)
	@echo " ---> Building $@ ... " 
ifeq ($(ARCH),aix)
		@/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		@/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
          $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif 
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "  ---> $@ done"

#
# clean
#
clean:
	@$(RM) *Dict.$(SrcSuf) *Dict.h *.o *.la *.lo *.dylib *.so *~

#
# generic 
#
.$(SrcSuf).$(ObjSuf):
	@echo "  ---> Compiling $@ ... "
	$(CXX) $(CXXFLAGS) -c $<
	@echo "  ---> $@ done"
#
# rules to get the ROOT dictionnary for each individual file
#
%Dict.C: %.h 
	@echo "  ---> Generating ROOT dictionary $@ ... "
	@rootcint -f $*Dict.C -c -I./ $*.h $*LinkDef.h
	@echo "  ---> $@ done"


