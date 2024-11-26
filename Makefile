# Set compiler and linker flags
CC			= g++
AR			= ar

# NOTE: If the dependencies are not installed in the default locations, adjust the path.
PREFIX		= $(HOME)
CFLAGS  	= -Wall -pthread -std=gnu++0x -Ofast
CPPFLAGS 	= -I $(PREFIX)/openssl-3.3.0/include
LNFLAGS  	= -Wl,-rpath=$(PREFIX)/openssl-3.3.0/lib64 -L $(PREFIX)/openssl-3.3.0/lib64 -lntl -lgmp -lcrypto

SRCS		= $(wildcard *.cc)
OBJS		= $(SRCS:.cc=.o)

.PHONY: 	all clean

all: 		BLNS

BLNS:		$(OBJS)
			$(CC) $(CFLAGS) -o BLNS $(OBJS) $(CPPFLAGS) $(LNFLAGS)

%.o: 		%.cc params.h
			$(CC) $(CFLAGS) $(CPPFLAGS) -c $< 

clean:
			rm -f *.o
			rm -f BLNS