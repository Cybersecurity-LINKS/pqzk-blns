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

# Enable or disable Falcon keygen and GSampler
USE_FALCON 	= 1
#USE_FALCON = 0

# Set compiler and linker flags
CC			= g++
CFLAGS  	= -Wall -pthread -std=gnu++0x -Ofast
LNFLAGS  	= -lntl -lgmp 

# If necessary, set Falcon library location and flags
ifeq ($(USE_FALCON),1)
FALCON_DIR 	= /Falcon-impl-20211101
FALCON_PATH	= $(CURDIR)$(FALCON_DIR)/
FALCON_LIB  = $(FALCON_PATH)/libfalcon.a
FILE_EXISTS = $(or $(and $(wildcard $(FALCON_LIB)),1),0)

CFLAGS  	+= -Drestrict=__restrict__ -DENABLE_FALCON -I$(FALCON_PATH)
LNFLAGS 	+= -L$(FALCON_PATH) -lfalcon
endif

SRCS		= $(wildcard *.cc)
OBJS		= $(SRCS:.cc=.o)


.PHONY: 	all falcon clean cleanall

all: 		falcon BLNS

falcon:	
ifeq ($(USE_FALCON),1)	
ifeq ($(FILE_EXISTS),0) # Falcon files do not already exists in $(FALCON_PATH)
	$(info Download and build Falcon...)
	wget https://falcon-sign.info/Falcon-impl-20211101.zip;\
	unzip Falcon-impl-20211101.zip;\
	cd Falcon-impl-20211101;\
	make -j$(nproc);\
	ar rcs libfalcon.a codec.o common.o falcon.o fft.o fpr.o keygen.o rng.o shake.o sign.o vrfy.o;\
	cd ..
endif
endif

BLNS:		$(OBJS)
			$(CC) $(CFLAGS) -o BLNS $(OBJS) $(LNFLAGS)

%.o: 		%.cc params.h
			$(CC) $(CFLAGS) -c $< 

clean:
			rm -f *.o
			rm -f BLNS

cleanall:   clean
ifeq ($(USE_FALCON),1)	
			rm -rf $(FALCON_PATH)
endif
