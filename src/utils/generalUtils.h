#ifndef GENERAL_UTILS_H
#define GENERAL_UTILS_H
/**
 * @file general.h
 * @author Caleb Hash
 * @date March 01, 2025
 * @brief General utilities library
 */

#include <chrono>

class Timer {
  private:
    std::chrono::high_resolution_clock::time_point start_time, end_time;
    bool running = false;

  public:
    void start();
    void end();
    void result(const char* message) const;
};
#endif // GENERAL_UTILS_H
