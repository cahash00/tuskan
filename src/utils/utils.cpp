/**
 * @file utils.cpp
 * @brief utilities library
 * @details will contain general utilities like timers, etc...
 * @author Caleb Hash
 * @date 2025-02-23
 */

/**
 * @brief Timer class
 */
class Timer {
  private:
    std::chrono::high_resolution_clock::time_point start_time, end_time;
    bool running = false;

  public:
    void start() {
      start_time = std::chrono::high_resolution_clock::now();
      running = true;
    }

    void end() {
      if (running) {
        end_time = std::chrono::high_resolution_clock::now();
        running = false;
      } else {
        std::cerr << "Timer was not started!\n";
      }
    }

    void result(const char* message) const {
      if (!running) {
        std::chrono::duration<double> elapsed = end_time - start_time;
        std::cout << " " << scientific << setprecision(15) <<message << " " << elapsed.count() << endl;
      } else {
        std::cerr << "Timer is still running! Call end() before result().\n";
      }
    }
};
