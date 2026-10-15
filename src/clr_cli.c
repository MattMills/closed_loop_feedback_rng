/*
 * clr -- stream generator output.
 *
 *   clr [--single] [--seed TEXT | --seed-hex HEX] [--rounds N] [--bytes N | --hex N]
 *
 * Default is the dual generator (clr_dual_*); --single uses one bare CLR loop.
 * With no --bytes/--hex, raw bytes are written to stdout until the pipe closes
 * (e.g. `clr | RNG_test stdin`).  With no seed, 32 bytes from /dev/urandom.
 */
#include "clr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void)
{
    fprintf(stderr, "usage: clr [--single] [--seed TEXT | --seed-hex HEX] [--rounds N] [--bytes N | --hex N]\n");
    exit(2);
}

static int hexval(int c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

int main(int argc, char **argv)
{
    uint8_t seed[4096];
    size_t seed_len = 0;
    int have_seed = 0, hex = 0, single = 0;
    unsigned long rounds = 0;
    unsigned long long count = 0;
    int bounded = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--single")) {
            single = 1;
        } else if (!strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed_len = strlen(argv[++i]);
            if (seed_len > sizeof seed) usage();
            memcpy(seed, argv[i], seed_len);
            have_seed = 1;
        } else if (!strcmp(argv[i], "--seed-hex") && i + 1 < argc) {
            const char *h = argv[++i];
            size_t n = strlen(h);
            if (n % 2 || n / 2 > sizeof seed) usage();
            for (size_t k = 0; k < n / 2; k++) {
                int a = hexval(h[2 * k]), b = hexval(h[2 * k + 1]);
                if (a < 0 || b < 0) usage();
                seed[k] = (uint8_t)(a << 4 | b);
            }
            seed_len = n / 2;
            have_seed = 1;
        } else if (!strcmp(argv[i], "--rounds") && i + 1 < argc) {
            rounds = strtoul(argv[++i], NULL, 10);
            if (rounds == 0) usage();
        } else if ((!strcmp(argv[i], "--bytes") || !strcmp(argv[i], "--hex")) && i + 1 < argc) {
            hex = !strcmp(argv[i], "--hex");
            count = strtoull(argv[++i], NULL, 10);
            bounded = 1;
        } else {
            usage();
        }
    }

    if (!have_seed) {
        FILE *f = fopen("/dev/urandom", "rb");
        if (!f || fread(seed, 1, 32, f) != 32) {
            fprintf(stderr, "clr: cannot read /dev/urandom\n");
            return 1;
        }
        fclose(f);
        seed_len = 32;
    }

    static clr_state st;
    static clr_dual dual;
    if (single)
        clr_init(&st, seed, seed_len, (uint32_t)rounds);
    else
        clr_dual_init(&dual, seed, seed_len, (uint32_t)rounds);
    memset(seed, 0, sizeof seed);

    uint8_t buf[1 << 16];
    if (hex) {
        while (count) {
            size_t n = count < sizeof buf ? (size_t)count : sizeof buf;
            if (single) clr_fill(&st, buf, n); else clr_dual_fill(&dual, buf, n);
            for (size_t k = 0; k < n; k++)
                printf("%02x", buf[k]);
            count -= n;
        }
        printf("\n");
        return 0;
    }
    for (;;) {
        size_t n = sizeof buf;
        if (bounded) {
            if (!count) break;
            if (count < n) n = (size_t)count;
            count -= n;
        }
        if (single) clr_fill(&st, buf, n); else clr_dual_fill(&dual, buf, n);
        if (fwrite(buf, 1, n, stdout) != n)
            return 0;                           /* reader closed the pipe */
    }
    return 0;
}
