#include "StarExecutable.hpp"
#include <iostream>
#include <string_view>

int main(int argc, char** argv) {
    const bool expectDisabled = argc == 2 && std::string_view(argv[1]) == "--expect-disabled";
    try {
        const auto path = argc == 2 && !expectDisabled
                              ? utility::star::executablePath(argv[1])
                              : utility::star::executablePath();
        std::cout << path.string() << '\n';
        return expectDisabled ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return expectDisabled && !utility::star::enabled() ? 0 : 1;
    }
}
