CC     ?= cc
CFLAGS ?= -O2 -Wall -Wextra -Wpedantic -std=c99
BUILD  := build
HDRS   := src/clr.h src/clr_tables.h
LIB    := src/clr.c src/clr_dual.c src/clr_entropy.c

all: $(BUILD)/clr $(BUILD)/test_clr $(BUILD)/libclr.so $(BUILD)/clr_linear_control $(BUILD)/libclr_linear_control.so

$(BUILD):
	mkdir -p $@

$(BUILD)/clr: src/clr_cli.c $(LIB) $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -o $@ src/clr_cli.c $(LIB)

# Same architecture with every nonlinear step made linear -- a control for analysis/ only.
$(BUILD)/clr_linear_control: src/clr_cli.c $(LIB) $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -DCLR_LINEAR_CONTROL -o $@ src/clr_cli.c $(LIB)

$(BUILD)/libclr_linear_control.so: $(LIB) $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -DCLR_LINEAR_CONTROL -fPIC -shared -o $@ $(LIB)

$(BUILD)/test_clr: tests/test_clr.c $(LIB) $(HDRS) tests/kat_abc.inc tests/kat_dual_abc.inc | $(BUILD)
	$(CC) $(CFLAGS) -o $@ tests/test_clr.c src/clr_entropy.c

$(BUILD)/libclr.so: $(LIB) $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -fPIC -shared -o $@ $(LIB)

test: $(BUILD)/test_clr
	./$(BUILD)/test_clr
	python3 tools/gen_tables.py --check

tables:
	python3 tools/gen_tables.py

clean:
	rm -rf $(BUILD)

.PHONY: all test tables clean
