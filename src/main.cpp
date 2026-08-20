// Standard
#include <execinfo.h>
#include <unistd.h>

#include <csignal>
#include <cstdlib>

// Boost
#include <boost/program_options.hpp>

// Internal
#include "AnnotationHierarchyError.hpp"
#include "Config.hpp"
#include "Logger.hpp"
#include "Runner.hpp"
#include "UnderlyingSequence.hpp"
#include "Utility.hpp"

auto main(int argc, const char* const argv[]) -> int {
    signal(SIGSEGV, helper::crashHandler);
    try {
        Runner::runPipeline(argc, argv);
        return EXIT_SUCCESS;
    } catch (const annotation::AnnotationHierarchyError& error) {
        Logger::log<SourceLocation{}, LogLevel::ERROR>(error.what());
    }
    return EXIT_FAILURE;
}
