### Compiler options
CC 		 = gcc
CPP 	 = g++
CFLAGS   = -O3 -DNDEBUG -Wall -g 
CPPFLAGS = -O3 -DNDEBUG -Wall -g -std=c++11

LDFLAGS  = 
LIBS     =  -lm
LINKFLAGS = -Wl,-rpath,libdir
includedir   = /usr/local/sundials-6/include
libdir       = /usr/local/sundials-6/lib


CPPSOURCES = odesolver 
CPPSOURCES_DEPENDENCIES = 
OBJECTS = ${CPPSOURCES:=.o}
OBJECTS_DEPENDENCIES = ${CPPSOURCES_DEPENDENCIES:=.o}
DEPS = ${CPPSOURCES:=.d}

OBJDIR = obj
DEPDIR = dep




##########################################################################################
##################### Shouldn't need to change from here downwards! ######################
##########################################################################################

# -------------- include headers and libaries from Sundials CVODE package --------------- #
TMP_INCS  = ${includedir} ${INCLUDES_SLUMT} ${INCLUDES_KLU}
INCLUDES  = $(addprefix -I, ${TMP_INCS})
LIBRARIES = -lsundials_cvode -lsundials_nvecserial -lsundials_nvecmanyvector ${LIBS}
# --------------------------------------------------------------------------------------- #


# ------------------------- rules to make object files ---------------------------------- #
${OBJDIR}/%.o: %.cpp
	@mkdir -p ${OBJDIR}
	${CPP} ${CPPFLAGS} ${CFLAGS} ${INCLUDES}  -o $@ -c $<

### rules to make dependency files
${DEPDIR}/%.d: %.cpp
	@mkdir -p ${DEPDIR}
	${CPP} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} $< -MM -MT $(@:%.d=%.o) >$@
# -------------------------------------------------------------------------------------- #


# ---------------------------- Targets beginning here ---------------------------------- #
all: targets 

targets: ${OBJDIR}/${OBJECTS} ${DEPDIR}/${DEPS}
	@for i in ${CPPSOURCES} ; do \
	  echo "${CPP} -o $${i} ${OBJDIR}/$${i}.o ${OBJECTS_DEPENDENCIES} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CPP} -o $${i} ${OBJDIR}/$${i}.o ${OBJECTS_DEPENDENCIES} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done



${OBJECTS}: ${OBJECTS_DEPENDENCIES}
### could create action here if wanted



# ---------------------------------------- Cleaning ------------------------------------ #

# Cleans complete project
.PHONY: clean`
clean:
	rm -f ${OBJECTS_DEPENDENCIES}
	rm -f ${OBJDIR}/${OBJECTS}
	rm -f ${DEPDIR}/${DEPS}
	rm -f ${CPPSOURCES}

# Cleans only all files with the extension .d
.PHONY: cleandep
cleandep:
	rm -f ${DEPDIR}/${DEPS}

# Cleans everything including removing csv datafiles, obj and dep directories
.PHONY: cleanALL
cleanALL:
	rm -f ${OBJECTS_DEPENDENCIES}
	rm -rf ${OBJDIR}
	rm -rf ${DEPDIR}
	rm -f ${CPPSOURCES}
	rm -f *_sol.csv
	rm -f *_err.csv
	rm -f *_stats.csv

# --------------------------------------------------------------------------------------- #