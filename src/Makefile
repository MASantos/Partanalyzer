#/*
#*    This file is part of Partanalyzer.
#*
#*    Partanalyzer is free software: you can redistribute it and/or modify
#*    it under the terms of the GNU General Public License as published by
#*    the Free Software Foundation, either version 3 of the License, or
#*    (at your option) any later version.
#*
#*    Partanalyzer is distributed in the hope that it will be useful,
#*    but WITHOUT ANY WARRANTY; without even the implied warranty of
#*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*    GNU General Public License for more details.
#*
#*    You should have received a copy of the GNU General Public License
#*    along with Partanalyzer.  If not, see <http://www.gnu.org/licenses/>.
#*
#*/
#/** Partanalyzer
#Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2009.
#Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
#*/

SHELL = /bin/sh

WDIR = ..
SRCD = src
srcdir = $(WDIR)/$(SRCD)

.SUFFIXES:
.SUFFIXES: .cc .cpp .o .h

CC=g++
CPPFLAGS =-O2
LDFLAGS = 
#CPPFLAGS="-DDEBUG"
BIN = partanalyzer
SRC = \
	Roulette.cc \
	BellNumber.cc \
	Sequence.cc \
	MultipleSeqAlign.cc \
	sNeighborhood.cc \
	Partition.cc \
	PartitionStats.cc \
	MatrixOfValues.cc \
	Ccop.cc \
	partanalyzer_help.cc \
	partanalyzer.cc 

OBJ = \
	Roulette.o \
	BellNumber.o \
	Sequence.o \
	MultipleSeqAlign.o \
	sNeighborhood.o \
	Partition.o \
	PartitionStats.o \
	MatrixOfValues.o \
	Ccop.o \
	partanalyzer_help.o \
	partanalyzer.o 

HEADERS = \
	Roulette.h \
 	BellNumber.h \
	Sequence.h \
	MultipleSeqAlign.h \
	sNeighborhood.h \
	Partition.h \
	PartitionStats.h \
	MatrixOfValues.h \
	Ccop.h \
	partanalyzer_help.h \
	partanalyzer_includes.h \
	partanalyzer_definitions.h \
	partanalyzer_basic_operations.h \
	partanalyzer.h 

SRCS = $(SRC) $(HEADERS)

TAGS=tags
CTAGS=ctags > $(TAGS)
DOC_BIN=doxygen
DOC_CONF=$(WDIR)/Doxyfile
DOC_DIR=$(WDIR)/Docs

CPTOUSERS = $(WDIR)/cp2users

TESTSDIR=$(WDIR)/Tests
TCONSENS=$(TESTSDIR)/Consensus
TCONFAIL=$(TCONSENS)/Fail
TCONOK=$(TCONSENS)/OK
THASSE=$(TESTSDIR)/Hasse

LINKBIN = [ -L $(WDIR)/$(BIN) ] || ln -s $(SRCD)/$(BIN) $(WDIR)/$(BIN)

#############################################   B U I L D    #####################################

default: $(TAGS) $(BIN)

copy:	default $(CPTOUSERS)
	$(CPTOUSERS)

touch:	
	touch $(SRCS)

debug: $(TAGS) touch 
	CPPFLAGS="$(CPPFLAGS) -DDEBUG" make 
$(TAGS):
	$(CTAGS) $(SRCS)

$(BIN): $(OBJ) $(HEADERS)
	$(CC) -o $@ $(LDFLAGS) $(OBJ)
	$(LINKBIN)

$(BIN).DEBUG: $(OBJ) $(HEADERS)
	$(CC) -o $(BIN) $(CPPFLAGS) -DDEBUG $(OBJ)
	$(LINKBIN)

compile: $(OBJ)
	
%.cc%.o%.h:  
	$(CC) -g $(CPPFLAGS) $<

clean: 
	rm $(OBJ) 2>/dev/null
#	[ `(ls *.o | wc -l ) 2>/dev/null` -gt 0 ] && rm $(BIN) $(OBJ) $(TAGS) 2>/dev/null

distclean: 
	[ `(ls *.o | wc -l ) 2>/dev/null` -gt 0 ] && rm $(BIN) $(OBJ) $(TAGS) 2>/dev/null

documentation:
	$(DOC_BIN) $(DOC_CONF)
	cd $(DOC_DIR)/latex  && make 

#############################################   B U I L D    #####################################

#############################################   T E S T S   #####################################
tests: $(BIN) start_tests_msg test_consensus test_hasse

start_tests_msg:
	
test_consensus: start_cons_msg test_cons_ok test_cons_fail

start_cons_msg:
	echo "Testing consensus"
test_cons_ok:
	cd $(TCONOK) && make -s BIN=$(BIN) BINDIR=$(WDIR) VERBOSE=$(TVERB)

test_cons_fail:
	cd $(TCONFAIL) && make -s BIN=$(BIN) BINDIR=$(WDIR) VERBOSE=$(TVERB)

start_hasse_msg:
test_hasse: start_hasse_msg
	cd $(THASSE) && make -s BIN=$(BIN) BINDIR=$(WDIR) VERBOSE=$(TVERB)

#############################################   T E S T S   #####################################
