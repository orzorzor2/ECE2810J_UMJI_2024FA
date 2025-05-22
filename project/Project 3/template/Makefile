CC = clang++
CFLAGS = -std=c++1z -O2 -Wall -Wextra -Werror -pedantic -Wno-unused-result -Wconversion -Wvla

TARGET = main

all: $(TARGET)

$(TARGET): main.cpp
	$(CC) $(CFLAGS) main.cpp -o $(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)
