
CC = gcc
CFLAGS = -O2 -Wall -I../libbayesrco_stats -I../libbayesrco_io
LDFLAGS = -L../libbayesrco_stats -L../libbayesrco_io -lbayesrco_io -lbayesrco_stats -lm

TARGET = clibayesrco
SOURCES = main.c
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
