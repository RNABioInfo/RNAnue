// Standard
#include <execinfo.h>
#include <unistd.h>

#include <csignal>
#include <cstdlib>

// Boost
#include <boost/program_options.hpp>

// Internal
#include "Config.hpp"
#include "Runner.hpp"
#include "Utility.hpp"

using namespace pipelines::analyze;

auto main(int argc, const char* const argv[]) -> int {
    signal(SIGSEGV, helper::crashHandler);

    Runner::runPipeline(argc, argv);
}
