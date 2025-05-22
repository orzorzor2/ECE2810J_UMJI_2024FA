CC = g++
CFLAGS = -std=c++17 -Wall -Wextra -Werror -pedantic -Wno-unused-result -Wconversion -Wvla
DEBUG_FLAGS = -g -DDEBUG

TARGET = main
SRC = main.cpp
HEADERS = sokoban.hpp
TEST_DIR = test_cases
OUTPUT_DIR = outputs

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

debug: CFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

test_all: $(TARGET)
	@for file in $(TEST_DIR)/*.in; do \
		echo "Testing $${file}..."; \
		./$(TARGET) < $${file} > $(OUTPUT_DIR)/$$(basename $${file} .in).out 2> $(OUTPUT_DIR)/$$(basename $${file} .in).err; \
	done

test_single: $(TARGET)
	@echo "Testing $(TEST_CASE)"
	./$(TARGET) < $(TEST_DIR)/$(TEST_CASE).in

memory_check: $(TARGET)
	valgrind --leak-check=full ./$(TARGET) < $(TEST_DIR)/$(TEST_CASE).in

code_check:
	cpplint --linelength=120 --filter=-legal,-readability/casting,-whitespace,-runtime/printf,-runtime/threadsafe_fn,-readability/todo,-build/include_subdir,-build/header_guard *.h *.cpp
	cppcheck --enable=all *.cpp

setup:
	@mkdir -p $(TEST_DIR)
	@mkdir -p $(OUTPUT_DIR)

clean:
	rm -f $(TARGET)
	rm -rf $(OUTPUT_DIR)/*

.PHONY: all debug test_all test_single memory_check code_check setup clean