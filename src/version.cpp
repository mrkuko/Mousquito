#include "version.hpp"

#ifndef MQ_VERSION
#define MQ_VERSION "0.0.0"
#endif

namespace mq {

const char* version() { return MQ_VERSION; }

} // namespace mq
