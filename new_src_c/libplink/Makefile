CC = gcc
CFLAGS = -O2 -Wall -fPIC
AR = ar
ARFLAGS = rcs

SOURCES = plink_io.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = libplink.a

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
