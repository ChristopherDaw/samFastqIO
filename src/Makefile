# Makefile for building C programs

SRC=main.c sam_models.c sam_file_allocation.c qv_codebook.c qv_compression.c read_compression.c read_decompression.c util.c well.c Arithmetic_stream.c stream_model.c pmf.c distortion.c quantizer.c sfio_ssh_session.c

OBJ=$(SRC:.c=.o)

CC=gcc
RM=rm -f

CFLAGS=-O3 -ffast-math -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -Wall -I../include -I/usr/local/include -L/usr/local/lib -lssh -D__APPLE__
LDFLAGS=-L/usr/local/lib -lc -lm -lssh

%.o : %.c
	$(CC) $(CFLAGS) -c $<

all : sfio

sfio : $(OBJ)
	$(CC) $(OBJ) -o sfio $(LDFLAGS)

debug : CFLAGS += -DDEBUG -ggdb
debug : sfio

clean :
	$(RM) *.o sfio
