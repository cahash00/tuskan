/**
 * @file general.cpp
 * @author Caleb Hash
 * @date March 01, 2025
 * @brief General utilities library
 */

#include <generalUtils.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <spdlog/spdlog.h>

void Timer::start() {
  start_time = std::chrono::high_resolution_clock::now();
  running = true;
}

void Timer::stop() {
  if (running) {
    end_time = std::chrono::high_resolution_clock::now();
    running = false;
  } else {
    std::cerr << "Timer was not started!\n";
  }
}

void Timer::result(const char* message) const {
  if (!running) {
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << std::scientific << std::setprecision(6) << message 
              << " " << elapsed.count() << " seconds" << std::endl;
  } else {
    std::cerr << "Timer is still running! Call end() before result().\n";
  }
}
double Timer::time() {
  std::chrono::duration<double> elapsed = end_time - start_time;
  double elapsedDouble = elapsed.count();
  return elapsedDouble;
}
/******************************************************************************/
/**
 * @brief Custom formatter for relative timestamps
 */
customSPDLOG::customSPDLOG() : start_time(std::chrono::steady_clock::now()) {}

void customSPDLOG::format(const spdlog::details::log_msg&,
                                   const std::tm&,
                                   spdlog::memory_buf_t& dest) {
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time);

    fmt::format_to(std::back_inserter(dest), "[{:03}.{:03}]", duration.count() / 1000, duration.count() % 1000);
}

std::unique_ptr<spdlog::custom_flag_formatter> customSPDLOG::clone() const {
    return std::make_unique<customSPDLOG>();
}

