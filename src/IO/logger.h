#ifndef LOGGER_H
#define LOGGER_H

#include <memory>
#include <generalUtils.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/pattern_formatter.h>

namespace IO {

extern std::shared_ptr<spdlog::logger> logger;

void init_logger();

} // end namespace IO
#endif // end LOGGER_H
