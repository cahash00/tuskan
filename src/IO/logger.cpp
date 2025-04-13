#ifndef LOGGER_H
#define LOGGER_H

#include <memory>
#include <generalUtils.h>
#include <logger.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/pattern_formatter.h>
#include <mpi.h>

namespace IO {

std::shared_ptr<spdlog::logger> logger = nullptr;

void init_logger() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // ... create the logger
  logger = spdlog::stdout_color_mt("rank_"+std::to_string(rank));
  
  // Apply a custom formatter with relative timestamps
  auto formatter = std::make_unique<spdlog::pattern_formatter>();
  formatter->add_flag<customSPDLOG>('R'); // Use %R for relative time
  formatter->set_pattern("%R [%^%l%$] %v"); // Example: [12.345] [INFO] Message
  logger->set_formatter(std::move(formatter));

  // Set the logger as default
  if (rank==0) {
    spdlog::set_default_logger(logger);
  }
} // end init_logger

void logit(const std::string& message) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank==0) {
    logger->info(message);
  }
}
  


} // end namespace IO
#endif // end LOGGER_H
