#ifndef LOGGER_H
#define LOGGER_H

#include <memory>
#include <generalUtils.h>
#include <logger.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/pattern_formatter.h>

namespace IO {

std::shared_ptr<spdlog::logger> logger = nullptr;

void init_logger() {
  // ... create the logger
  logger = spdlog::stdout_color_mt("console");
  // Apply a custom formatter with relative timestamps
  auto formatter = std::make_unique<spdlog::pattern_formatter>();
  formatter->add_flag<customSPDLOG>('R'); // Use %R for relative time
  formatter->set_pattern("%R [%^%l%$] %v"); // Example: [12.345] [INFO] Message
  logger->set_formatter(std::move(formatter));
  // Set the logger as default
  spdlog::set_default_logger(logger);
} // end init_logger
} // end namespace IO
#endif // end LOGGER_H
