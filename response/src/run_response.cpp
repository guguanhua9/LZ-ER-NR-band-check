#include "nestDiamxCore.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

// The only editable job control is population, event count, and file paths.
// The selected physics parameters, gain configuration, and seeds are frozen.
int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            std::cerr << "Usage: lz_response er|nr EVENTS SPECTRUM.csv OUTPUT.csv\n";
            return 2;
        }
        const std::string population = argv[1];
        if (population != "er" && population != "nr")
            throw std::runtime_error("population must be er or nr");
        const std::string count_text = argv[2];
        if (count_text.empty() || count_text.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("EVENTS must be a positive integer");
        const auto events = std::stoull(count_text);
        if (events == 0) throw std::runtime_error("EVENTS must be positive");
        if (!std::filesystem::is_regular_file(argv[3]))
            throw std::runtime_error("SPECTRUM.csv must exist as a regular file");
        // Open the output before starting the job, so a path error is immediate.
        std::ofstream output(argv[4]);
        if (!output) throw std::runtime_error("Unable to open OUTPUT.csv");
        const std::string interaction = population == "er" ? "beta" : "NR";
        const uint64_t seed = population == "er" ? 2452026 : 260903;
        const auto result = LZModel("lz_ws2024", interaction, seed, argv[3],
                                    events, 1.0, 1.0, 0);
        if (result.s1c_phd.size() != events || result.s2c_phd.size() != events ||
            result.energy_rec.size() != events)
            throw std::runtime_error("Unexpected output vector lengths");
        output << std::setprecision(std::numeric_limits<double>::max_digits10)
               << "s1c_phd,s2c_phd,energy_rec\n";
        for (size_t i = 0; i < events; ++i)
            output << result.s1c_phd[i] << ',' << result.s2c_phd[i] << ','
                   << result.energy_rec[i] << '\n';
        if (!output) throw std::runtime_error("Output write failed");
        std::cerr << "Wrote " << events << " events for " << interaction
                  << "; seed=" << seed << "; er_omega_scale=1; ablation=0\n";
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "lz_response: " << exc.what() << '\n';
        return 1;
    }
}
