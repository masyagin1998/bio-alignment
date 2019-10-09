SRC_PREFIX=src/
BIN_PREFIX=bin/
OBJS_PREFIX=.objs/

CC=gcc

CFLAGS=-g -Wall -Wextra -std=c89 -O3

all: $(BIN_PREFIX)nw

FASTA_SRC_PREFIX=$(SRC_PREFIX)fasta/
FASTA_SRC=$(shell find $(FASTA_SRC_PREFIX) -maxdepth 1 -name '*.c')
FASTA_OBJS_PREFIX=$(OBJS_PREFIX)fasta/
FASTA_OBJS=$(patsubst $(FASTA_SRC_PREFIX)%.c,$(FASTA_OBJS_PREFIX)%.o,$(FASTA_SRC))

$(FASTA_OBJS_PREFIX)%.o: $(FASTA_SRC_PREFIX)%.c $(FASTA_SRC_PREFIX)%.h
	mkdir -p $(FASTA_OBJS_PREFIX)
	$(CC) $(CFLAGS) -c $< -o $@

NW_SRC_PREFIX=$(SRC_PREFIX)needleman-wunsch/
NW_SRC=$(shell find $(NW_SRC_PREFIX) -maxdepth 1 -name '*.c')
NW_OBJS_PREFIX=$(OBJS_PREFIX)needleman-wunsch/
NW_OBJS=$(patsubst $(NW_SRC_PREFIX)%.c,$(NW_OBJS_PREFIX)%.o,$(NW_SRC))

$(NW_OBJS_PREFIX)%.o: $(NW_SRC_PREFIX)%.c $(NW_SRC_PREFIX)needleman-wunsch.h
	mkdir -p $(NW_OBJS_PREFIX)
	$(CC) $(CFLAGS) -c $< -o $@

MAIN_SRC_PREFIX=$(SRC_PREFIX)
MAIN_SRC=$(shell find $(MAIN_SRC_PREFIX) -maxdepth 1 -name '*.c')
MAIN_OBJS_PREFIX=$(OBJS_PREFIX)
MAIN_OBJS=$(patsubst $(MAIN_SRC_PREFIX)%.c,$(MAIN_OBJS_PREFIX)%.o,$(MAIN_SRC))

$(MAIN_OBJS_PREFIX)%.o: $(MAIN_SRC_PREFIX)%.c
	mkdir -p $(MAIN_OBJS_PREFIX)
	$(CC) $(CFLAGS) -I$(FASTA_SRC_PREFIX) -I$(NW_SRC_PREFIX) -c $< -o $@

$(BIN_PREFIX)nw: $(FASTA_OBJS) $(NW_OBJS) $(MAIN_OBJS)
	mkdir -p $(BIN_PREFIX)
	$(CC) $(CFLAGS) $^ -o $@


.PHONY: clean

clean:
	rm -rf $(BIN_PREFIX) $(OBJS_PREFIX)
