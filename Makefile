CC     ?= cc
CFLAGS ?= -O2 -Wall -Wextra -Wpedantic -std=c99
BUILD  := build
HDRS   := src/clr.h src/clr_tables.h

all: $(BUILD)/clr $(BUILD)/test_clr $(BUILD)/libclr.so $(BUILD)/clr_linear_control $(BUILD)/libclr_linear_control.so

$(BUILD):
	mkdir -p $@

$(BUILD)/clr: src/clr_cli.c src/clr.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -o $@ src/clr_cli.c src/clr.c

# Same architecture with every nonlinear step made linear -- a control for analysis/ only.
$(BUILD)/clr_linear_control: src/clr_cli.c src/clr.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -DCLR_LINEAR_CONTROL -o $@ src/clr_cli.c src/clr.c

$(BUILD)/libclr_linear_control.so: src/clr.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -DCLR_LINEAR_CONTROL -fPIC -shared -o $@ src/clr.c

$(BUILD)/test_clr: tests/test_clr.c src/clr.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -o $@ tests/test_clr.c

$(BUILD)/libclr.so: src/clr.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -fPIC -shared -o $@ src/clr.c

test: $(BUILD)/test_clr
	./$(BUILD)/test_clr
	python3 tools/gen_tables.py --check

tables:
	python3 tools/gen_tables.py

clean:
	rm -rf $(BUILD)

.PHONY: all test tables clean
