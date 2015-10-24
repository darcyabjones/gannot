C = $(shell pwd)
DATA = $(C)/data
GENOME_FILE = MNH120.fasta
ISOLATE = $(word 1, $(subst ., ,$(GENOME_FILE)))

THREADS = 1
PYTHON = python
RIPCRAWL = bin/ripcrawl.py -t -n 0 -w 500 -i 20 $(1) > $(2)
W2B = bin/windows2bedtrack.py -n $(1) -i $(2) -o $(3)
TRF=trf-4.0.7b $(1) 2 7 7 80 10 50 2000 -d -h -ngs > $(2)
TRF2GFF=bin/trf2gff.py
legacy_blast=/home/darcy/bin/blast-2.2.26/bin
TPSI=/usr/local/transposonpsi/08222010/transposonPSI.pl $(1) nuc
NN2GFF=bin/nn2gff.py -m 1 $(1) > $(2)

NSEQS = $(shell grep -c '>' $(1))

## Define filenames and directories

DATA = data

RIPCRAWL_DIR = RIPCrawl
RIPCRAWL_EXTS = .ripcrawl.csv .ripcrawl.CRI.bed .ripcrawl.GC.bed .ripcrawl.CRI.track.bed .ripcrawl.GC.track.bed
RIPCRAWL_FILES = $(foreach e, $(RIPCRAWL_EXTS), $(addprefix $(RIPCRAWL_DIR)/, $(addsuffix $(e), $(notdir $(basename $(GENOME_FILE))))))

TRF_DIR = trf
TRF_EXTS = .trf.dat .trf.gff3
TRF_FILES = $(foreach e, $(TRF_EXTS), $(addprefix $(TRF_DIR)/, $(addsuffix $(e), $(notdir $(basename $(GENOME_FILE))))))

TPSI_DIR = tpsi
TPSI_EXTS = .TPSI.allHits .TPSI.allHits.chains .TPSI.allHits.chains.gff3 \
	.TPSI.allHits.chains.bestPerLocus .TPSI.allHits.chains.bestPerLocus.gff3
TPSI_FILES = $(foreach e, $(TPSI_EXTS), $(addprefix $(TPSI_DIR)/, $(addsuffix $(e), $(notdir $(basename $(GENOME_FILE))))))

NNBLOCK_DIR = NNBlocks
NNBLOCK_FILES = $(addprefix $(NNBLOCK_DIR)/, $(addsuffix .NNBlocks.gff3, $(notdir $(basename $(GENOME_FILE)))))


## Commands

all: ripcrawl trf tpsi nnblocks

ripcrawl: $(RIPCRAWL_FILES)
trf: $(TRF_FILES)
tpsi: $(TPSI_FILES)
nnblocks: $(NNBLOCK_FILES)


$(RIPCRAWL_DIR)/%.ripcrawl.csv: $(DATA)/%$(suffix $(GENOME_FILE))
		mkdir -p $(dir $@)
		$(call RIPCRAWL, $<, $@)
$(RIPCRAWL_DIR)/%.ripcrawl.CRI.bed: $(RIPCRAWL_DIR)/%.ripcrawl.csv
		mkdir -p $(dir $@)
		tail -n +2 $< | awk '$$9!="fish" {print $$1 "\t" $$2 "\t" $$3 "\tCRI\t" $$8}' > $@
$(RIPCRAWL_DIR)/%.ripcrawl.GC.bed: $(RIPCRAWL_DIR)/%.ripcrawl.csv
		mkdir -p $(dir $@)
		tail -n +2 $< | awk '$$9!="fish" {print $$1 "\t" $$2 "\t" $$3 "\tGC\t" $$9}' > $@
$(RIPCRAWL_DIR)/%.ripcrawl.GC.track.bed: $(RIPCRAWL_DIR)/%.ripcrawl.GC.bed
		mkdir -p $(dir $@)
		$(call W2B, GC, $<, $@)
$(RIPCRAWL_DIR)/%.ripcrawl.CRI.track.bed: $(RIPCRAWL_DIR)/%.ripcrawl.CRI.bed
		mkdir -p $(dir $@)
		$(call W2B, CRI, $<, $@)

$(TRF_DIR)/%.trf.dat: $(DATA)/%$(suffix $(GENOME_FILE))
	mkdir -p $(TRF_DIR)
	$(call TRF, $<, $@)
$(TRF_DIR)/%.trf.gff3: $(TRF_DIR)/%.trf.dat
	mkdir -p $(TRF_DIR)
	$(call TRF2GFF, $<, $@)

$(foreach e, $(TPSI_EXTS), $(addprefix $(TPSI_DIR)/, $(addsuffix $(e), %))): $(DATA)/%$(suffix $(GENOME_FILE))
	mkdir -p $(dir $@)
	export PATH=$(legacy_blast):$(PATH); cd $(dir $@); $(call TPSI, $<)

$(NNBLOCKS_DIR)/%.NNBlocks.gff3: $(DATA)/%$(suffix $(GENOME_FILE))
	mkdir -p $(NNBLOCKS_DIR)
	$(call NN2GFF, $<, $@)
