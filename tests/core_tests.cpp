#include <cstdio>
#include <cstring>

#include "version.hpp"

static int failures = 0;
#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

static void testVersion() {
    const char* v = mq::version();
    CHECK(v != nullptr);
    CHECK(std::strlen(v) >= 5); // "X.Y.Z"
}

int main() {
    testVersion();
    if (failures == 0) std::puts("all tests passed");
    return failures == 0 ? 0 : 1;
}
