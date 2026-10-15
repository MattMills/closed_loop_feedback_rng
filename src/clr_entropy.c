/* OS entropy for the occluded twin.  getentropy() is in glibc >= 2.25, macOS and the BSDs. */
#define _DEFAULT_SOURCE

#include "clr.h"

#include <stdlib.h>
#include <unistd.h>

void clr_os_entropy(void *buf, size_t len)
{
    uint8_t *p = buf;
    while (len) {
        size_t n = len < 256 ? len : 256;
        if (getentropy(p, n) != 0)
            abort();                            /* never continue without entropy */
        p += n;
        len -= n;
    }
}
