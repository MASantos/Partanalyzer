
EXE=$(BINDIR)/$(BIN)


INPARTS=p1 p2 p3 p4 p5
CBASENAME=consensus

default: $(INPARTS) startmsg genconsensus isitaPart
	printf "\t\t\t\t\t\t\t\t\t\tOK\n"

startmsg:
	printf "\tTesting valid consensus..."

genconsensus:
	$(EXE) $(VERBOSE) -clstat -cons $(INPARTS) > $(CBASENAME).out
	sed -n '/^#BeginConsensusPartition/,/^#EndConsensusPartition/p' $(CBASENAME).out | grep -v "^#" > $(CBASENAME)
	test "x`grep ERROR $(CBASENAME)|wc -l`" == "x0" || less $(CBASENAME)

isitaPart:
	$(EXE) $(VERBOSE) --is-partition $(CBASENAME) >  isPart-$(CBASENAME).out
	test "x`grep ERROR $(CBASENAME).out | wc -l`" == "x0" || less $(CBASENAME).out
