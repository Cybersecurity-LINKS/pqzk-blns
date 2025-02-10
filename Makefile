# Copyright 2025 Fondazione LINKS

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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