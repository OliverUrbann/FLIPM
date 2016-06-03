CC        = gcc
OS        = $(shell uname)
CFLAGS    = -Wall -I. -o boundedness
LDFLAGS   = -lm

ifneq "$(OS)" "AmigaOS"
CFLAGS    += -std=c99
endif

boundedness:	main.c boundedness.c boundedness.h
	gcc $(CFLAGS) main.c boundedness.c $(LDFLAGS)
