/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   #
#    Conflict Resolution for Assured Contingency Landing Management                 #
#    Copyright (C) 2026 Huseyin Emre Tekaslan                                       #
#                                                                                   #
#    This program is free software: you can redistribute it and/or modify           #
#    it under the terms of the GNU General Public License as published by           #
#    the Free Software Foundation, either version 3 of the License, or              #
#    (at your option) any later version.                                            #
#                                                                                   #
#    This program is distributed in the hope that it will be useful,                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  #
#    GNU General Public License for more details.                                   #
#                                                                                   #
#    You should have received a copy of the GNU General Public License              #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.         #
#                                                                                   #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Conflict Resolution for Assured Contingency Landing Management                 %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : July 2026                                                         %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
/**
 * @file main.cpp
 * @brief Generate conflict-resolution advisories for available planner outputs.
 *
 * Solution suffixes:
 *   0: Search
 *   1: Dubins
 *
 * A solution is processed only when both its emergency trajectory and
 * corresponding conflict log exist. If both Search and Dubins solutions
 * exist, both are processed.
 */

#include "conflict.hpp"

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

/**
 * Number of resolution exit-flag categories returned by Conflict::resolve().
 */
constexpr std::size_t kExitFlagCount = 8;

/**
 * Identifier passed to functions that still require a case number.
 *
 * The current PAA workflow contains one planning scenario rather than
 * numbered cases.
 */
constexpr std::size_t kRunIdentifier = 0;

/**
 * Description of one emergency-landing solution type.
 */
struct SolutionDefinition {
    PathType pathType;
    std::string name;
    fs::path emergencyPathFilename;
    fs::path conflictLogFilename;
};

/**
 * Planner outputs that may be processed.
 */
const std::array<SolutionDefinition, 2> kSolutions{{
    {
        PathType::Search,
        "Search",
        "airspaceRisk0.csv",
        "conflicts0.csv"
    },
    {
        PathType::Dubins,
        "Dubins",
        "airspaceRisk1.csv",
        "conflicts1.csv"
    }
}};

/**
 * @brief Configure the conflict-resolution solver.
 *
 * @param conflict Conflict object to configure.
 */
void configureConflictSolver(Conflict& conflict)
{
    readHeatMap(&conflict.prohibited);

    conflict.dubinsopt.trajopt.verbose = 0;
    conflict.dubinsopt.trajopt.geoopt.verbose = 0;
    conflict.dubinsopt.trajopt.geoopt.model = WGS84;

    conflict.opt = Conflict::Optimizer::global;
    conflict.population = 4;
    conflict.tlimit = 60.0;

    conflict.windSpeed = 0.0;
    conflict.windDirection = 0.0;
}


/**
 * @brief Return the trajectory path for a solution.
 *
 * @param plannerOutputDirectory Planner output directory.
 * @param solution Solution definition.
 * @return Full path to the emergency trajectory.
 */
fs::path emergencyPath(
    const fs::path& plannerOutputDirectory,
    const SolutionDefinition& solution)
{
    return plannerOutputDirectory /
           solution.emergencyPathFilename;
}

/**
 * @brief Return the conflict-log path for a solution.
 *
 * @param plannerOutputDirectory Planner output directory.
 * @param solution Solution definition.
 * @return Full path to the conflict log.
 */
fs::path conflictLogPath(
    const fs::path& plannerOutputDirectory,
    const SolutionDefinition& solution)
{
    return plannerOutputDirectory /
           solution.conflictLogFilename;
}

/**
 * @brief Check whether both required files exist for a solution.
 *
 * @param plannerOutputDirectory Planner output directory.
 * @param solution Solution definition.
 * @return True when both the trajectory and conflict log exist.
 */
bool solutionExists(
    const fs::path& plannerOutputDirectory,
    const SolutionDefinition& solution)
{
    return fs::is_regular_file(
               emergencyPath(plannerOutputDirectory, solution)) &&
           fs::is_regular_file(
               conflictLogPath(plannerOutputDirectory, solution));
}

/**
 * @brief Report a solution for which only one required file exists.
 *
 * A complete solution requires both:
 *  - airspaceRisk%d.csv
 *  - conflicts%d.csv
 *
 * @param plannerOutputDirectory Planner output directory.
 * @param solution Solution definition.
 */
void reportIncompleteSolution(
    const fs::path& plannerOutputDirectory,
    const SolutionDefinition& solution)
{
    const fs::path trajectory =
        emergencyPath(plannerOutputDirectory, solution);

    const fs::path conflictLog =
        conflictLogPath(plannerOutputDirectory, solution);

    const bool hasTrajectory =
        fs::is_regular_file(trajectory);

    const bool hasConflictLog =
        fs::is_regular_file(conflictLog);

    /*
     * Nothing needs to be reported if both files exist or both files are
     * absent.
     */
    if (hasTrajectory == hasConflictLog) {
        return;
    }

    std::cerr
        << "Warning: incomplete "
        << solution.name
        << " solution";

    if (!hasTrajectory) {
        std::cerr
            << "; missing "
            << trajectory.string();
    }

    if (!hasConflictLog) {
        std::cerr
            << "; missing "
            << conflictLog.string();
    }

    std::cerr << '\n';
}

/**
 * @brief Add returned resolution flags to the aggregate counts.
 *
 * @tparam ExitFlagContainer Container returned by Conflict::resolve().
 * @param exitFlags Resolution exit flags.
 * @param flagCounts Aggregate exit-flag counts.
 * @param solutionName Name of the solution being processed.
 */
template <typename ExitFlagContainer>
void updateFlagCounts(
    const ExitFlagContainer& exitFlags,
    std::array<std::size_t, kExitFlagCount>& flagCounts,
    const std::string& solutionName)
{
    for (const auto rawExitFlag : exitFlags) {
        const auto exitFlag =
            static_cast<std::size_t>(rawExitFlag);

        if (exitFlag >= flagCounts.size()) {
            std::cerr
                << "Warning: invalid exit flag "
                << rawExitFlag
                << " returned for "
                << solutionName
                << '\n';

            continue;
        }

        ++flagCounts[exitFlag];
    }
}

/**
 * @brief Process one available emergency-landing solution.
 *
 * A separate Conflict object is created for each solution to prevent
 * trajectory, conflict-log, or optimizer state from carrying over from the
 * Search solution to the Dubins solution or vice versa.
 *
 * @param plannerOutputDirectory Planner output directory.
 * @param solutionOutputDirectory Output directory for this solution.
 * @param solution Solution definition.
 * @param caseLog Case-level CSV log.
 * @param resolutionLog Resolution-level CSV log.
 * @param flagCounts Aggregate exit-flag counts.
 */
void processSolution(
    const fs::path& plannerOutputDirectory,
    const fs::path& solutionOutputDirectory,
    const SolutionDefinition& solution,
    std::ofstream& caseLog,
    std::ofstream& resolutionLog,
    std::array<std::size_t, kExitFlagCount>& flagCounts)
{
    std::cout
        << solution.name
        << " solution\n";

    fs::create_directories(solutionOutputDirectory);

    Conflict conflict;
    configureConflictSolver(conflict);

    conflict.loadEmergencyPath(
        emergencyPath(plannerOutputDirectory, solution));

    conflict.loadConflictLog(
        conflictLogPath(plannerOutputDirectory, solution));

    if (conflict.conflictCount() == 0) {
        constexpr std::size_t noConflictFlag = 6;
        ++flagCounts[noConflictFlag];

        RiskMetrics risk{};

        printCaseRow(
            caseLog,
            kRunIdentifier,
            solution.name + "-NoConflict",
            0,
            0,
            0,
            0,
            risk);

        std::cout
            << "\tNo conflict\n";

        return;
    }

    /*
     * kRunIdentifier is retained because Conflict::resolve() currently
     * requires a numeric case identifier. The current workflow has no
     * numbered cases, so zero is used.
     */
    const auto exitFlags = conflict.resolve(
        kRunIdentifier,
        solutionOutputDirectory,
        caseLog,
        resolutionLog);

    updateFlagCounts(
        exitFlags,
        flagCounts,
        solution.name);
}

/**
 * @brief Print aggregate conflict-resolution statistics.
 *
 * When both Search and Dubins trajectories exist, the summary includes
 * advisories generated for both trajectories.
 *
 * @param flags Number of occurrences of each exit flag.
 */
void printResolutionSummary(
    const std::array<std::size_t, kExitFlagCount>& flags)
{
    std::cout
        << '\n'
        << "Conflict-resolution summary\n"
        << "---------------------------\n";

    std::printf(
        "%-22s %zu\n",
        "Take-off halt:",
        flags[0]);

    std::printf(
        "%-22s %zu\n",
        "Speed regulation:",
        flags[1]);

    std::printf(
        "%-22s %zu\n",
        "Altitude regulation:",
        flags[2]);

    std::printf(
        "%-22s %zu\n",
        "Lateral diversion:",
        flags[3]);

    std::printf(
        "%-22s %zu\n",
        "Go around and hold:",
        flags[4]);

    std::printf(
        "%-22s %zu\n",
        "No solution found:",
        flags[5]);

    std::printf(
        "%-22s %zu\n",
        "No conflict:",
        flags[6]);

    std::printf(
        "%-22s %zu\n",
        "Init. in conflict:",
        flags[7]);
}

}  // namespace

int main()
{
    try {
        /*
         * The executable is expected to run from PAA/avoid.
         *
         * current_path()               -> PAA/avoid
         * current_path().parent_path() -> PAA
         */
        const fs::path avoidDirectory =
            fs::current_path();

        const fs::path paaDirectory =
            avoidDirectory.parent_path();

        const fs::path plannerOutputDirectory =
            paaDirectory /
            "plan" /
            "out";

        const fs::path advisoryDirectory =
            avoidDirectory /
            "advisory";

        const fs::path caseLogPath =
            advisoryDirectory /
            "00_cases.csv";

        const fs::path resolutionLogPath =
            advisoryDirectory /
            "00_resolutions.csv";

        if (!fs::is_directory(plannerOutputDirectory)) {
            throw std::runtime_error(
                "Planner output directory does not exist: " +
                plannerOutputDirectory.string());
        }

        fs::create_directories(advisoryDirectory);

        std::ofstream caseLog(caseLogPath);
        std::ofstream resolutionLog(resolutionLogPath);

        if (!caseLog.is_open()) {
            throw std::runtime_error(
                "Unable to open case log: " +
                caseLogPath.string());
        }

        if (!resolutionLog.is_open()) {
            throw std::runtime_error(
                "Unable to open resolution log: " +
                resolutionLogPath.string());
        }

        printCasesHeader(caseLog);
        printResolutionsHeader(resolutionLog);

        std::array<std::size_t, kExitFlagCount> flagCounts{};

        std::size_t processedSolutionCount = 0;

        for (const auto& solution : kSolutions) {
            reportIncompleteSolution(
                plannerOutputDirectory,
                solution);

            if (!solutionExists(
                    plannerOutputDirectory,
                    solution)) {

                std::cout
                    << solution.name
                    << " solution not available\n";

                continue;
            }

            const fs::path solutionOutputDirectory =
                advisoryDirectory /
                solution.name;

            processSolution(
                plannerOutputDirectory,
                solutionOutputDirectory,
                solution,
                caseLog,
                resolutionLog,
                flagCounts);

            ++processedSolutionCount;
        }

        if (processedSolutionCount == 0) {
            throw std::runtime_error(
                "No complete Search or Dubins solution was found in: " +
                plannerOutputDirectory.string());
        }

        printResolutionSummary(flagCounts);

        std::cout
            << '\n'
            << "Solutions processed: "
            << processedSolutionCount
            << '\n'
            << "Advisory output: "
            << advisoryDirectory.string()
            << '\n';

        return EXIT_SUCCESS;
    }
    catch (const std::exception& exception) {
        std::cerr
            << "Error: "
            << exception.what()
            << '\n';

        return EXIT_FAILURE;
    }
}