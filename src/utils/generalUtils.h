#ifndef GENERAL_UTILS_H
#define GENERAL_UTILS_H
/**
 * @file general.h
 * @author Caleb Hash
 * @date March 01, 2025
 * @brief General utilities library
 */

#include <types.h>
#include <chrono>
#include <spdlog/spdlog.h>
#include <matar.h>

using namespace std;
using namespace mtr;
/**
 * @brief Timer class for timing things
 */
class Timer {
  private:
    chrono::high_resolution_clock::time_point start_time, end_time;
    bool running = false;

  public:
    void start();
    void stop();
    void result(const char* message) const;
    double time();
};
/******************************************************************************/
/**
 * @brief Custom formatter for relative timestamps
 */
class customSPDLOG : public spdlog::custom_flag_formatter {
   private:
     std::chrono::steady_clock::time_point start_time;

   public:
     customSPDLOG();

     void format(const spdlog::details::log_msg&,
         const std::tm&,
         spdlog::memory_buf_t& dest) override;

     std::unique_ptr<custom_flag_formatter> clone() const override;
 };
void print_field(FMatrix<double>& q);
#endif // GENERAL_UTILS_H
