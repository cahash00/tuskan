#ifndef ARGS_H
#define ARGS_H

#include <argparse/argparse.hpp>

namespace IO {
void getUserInput(int argc, char* argv[], argparse::ArgumentParser& program);
} // end namespace IO


#endif // end ARGS_H
