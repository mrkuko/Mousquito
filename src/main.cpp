#include <cstdio>
#include <cstring>

#include "version.hpp"

namespace {

void printUsage() {
    std::printf(
        "Mousquito %s - human-like mouse movement\n"
        "\n"
        "usage: mousquito <x,y> [x,y ...] [options]\n"
        "\n"
        "Not implemented yet: this build only prints this text.\n"
        "Hold Esc to swat it.\n",
        mq::version());
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "--version") == 0) {
        std::printf("mousquito %s\n", mq::version());
        return 0;
    }
    printUsage();
    return 0;
}
