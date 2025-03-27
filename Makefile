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
CFLAGS  	= -Wall -pthread -std=gnu++0x -Ofast -Drestrict=__restrict__
LNFLAGS  	= -lntl -lgmp -lfalcon  # Link Falcon library

# Falcon library location
FALCON_DIR	= $(HOME)/Falcon-impl-20211101/
FALCON_LIB  = $(FALCON_DIR)/libfalcon.a
FALCON_INC  = $(HOME)/Falcon-impl-20211101/  # Adjust if Falcon headers are in a different path


SRCS		= $(wildcard *.cc)
OBJS		= $(SRCS:.cc=.o)

.PHONY: all clean

all: BLNS

BLNS: $(OBJS) $(FALCON_LIB)
	$(CC) $(CFLAGS) -I$(FALCON_INC) -o BLNS $(OBJS) -L$(FALCON_DIR) $(LNFLAGS)

%.o: %.cc params.h
	$(CC) $(CFLAGS) -I$(FALCON_INC) -c $< 

clean:
	rm -f *.o
	rm -f BLNS
