#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "input.hpp"

using std::cerr;
using std::endl;

namespace fs = std::filesystem;

void parse_executable_arguments(
        PC_tree_t& conf_voicexx,
        long int& iter_start,
        int argc,
        char** argv,
        char const* const params_yaml)
{
    iter_start = 0;
    if (argc == 2) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            file.close();
            std::exit(EXIT_SUCCESS);
        }
    } else if (argc == 4) {
        if (argv[1] == std::string_view("--iter-restart")) {
            iter_start = std::strtol(argv[2], NULL, 10);
            conf_voicexx = PC_parse_path(fs::path(argv[3]).c_str());
        }
    } else {
        cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << endl;
        cerr << "or to perform a restart" << argv[0] << " [--iter-restart] <iter> <config_file.yml>"
             << endl;
        std::exit(EXIT_FAILURE);
    }
}

PC_tree_t parse_executable_arguments(int argc, char** argv, char const* const params_yaml)
{
    if (argc == 2) {
        return PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            file.close();
            std::exit(EXIT_SUCCESS);
        }
    }

    cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << endl;
    std::exit(EXIT_FAILURE);
}
