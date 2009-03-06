
EXE=$(BINDIR)/$(BIN)


INPARTS=p0  p1  p21  p22  p23  p24  p25  p26  p27  p31  p32  p33  p34  p35  p36
CBASENAME=hasse

default: $(INPARTS) startmsg getHasse
	printf "\t\t\t\t\t\t\t\t\t\t\tOK\n"

startmsg:
	printf "Testing hasse diagram..."

getHasse: gethasse_full

gethasse_full:
	$(EXE) $(VERBOSE) -hasse $(INPARTS) > $(CBASENAME).out
	sed -n '/^#BeginHasseNodes/,/^#EndHasseNodes/p' $(CBASENAME).out > $(CBASENAME).hierachy
	sed -n '/^#BeginHasseDiagram/,/^#EndHasseDiagram/p' $(CBASENAME).out > $(CBASENAME).graph

