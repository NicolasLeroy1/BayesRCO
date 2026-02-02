CC = gcc
CFLAGS = -O2 -Wall
LDFLAGS = -L../libbayesrco -L../libplink -lbayesrco -lplink -lm

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
