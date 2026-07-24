#include "aclm.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

/*
 * Metrics written to the result log for either the search-based path
 * or the minimum-risk Dubins path.
 */
typedef struct {
    double groundRisk;
    double airspaceRisk;
    double dGamma;
    double totalRuntime;
    double riskRuntime;
} PathMetrics;


/*
 * Create the output directory if it does not already exist.
 *
 * Returns:
 *   0 on success
 *  -1 on failure
 */
static int ensureDirectoryExists(const char *directory) {
    struct stat directoryInfo;

    if (stat(directory, &directoryInfo) == 0) {
        if (S_ISDIR(directoryInfo.st_mode)) {
            return 0;
        }

        fprintf(
            stderr,
            "Output path exists but is not a directory: %s\n",
            directory);

        return -1;
    }

    if (errno != ENOENT) {
        perror("Unable to inspect output directory");
        return -1;
    }

    if (mkdir(directory, 0775) != 0) {
        perror("Unable to create output directory");
        return -1;
    }

    return 0;
}


/*
 * Return unavailable path metrics for cases in which a path cannot
 * be generated or evaluated.
 */
static PathMetrics unavailableMetrics(void) {
    return (PathMetrics) {
        .groundRisk   = NAN,
        .airspaceRisk = NAN,
        .dGamma       = NAN,
        .totalRuntime = NAN,
        .riskRuntime  = NAN
    };
}


/*
 * Append one test-case result to the CSV log.
 *
 * Column order:
 *   1.  Exit flag
 *   2.  Initial latitude
 *   3.  Initial longitude
 *   4.  Initial altitude
 *   5.  Initial heading
 *   6.  Search ground risk
 *   7.  Search airspace risk
 *   8.  Search flight-path-angle change
 *   9.  Search runtime
 *   10. Search risk-evaluation runtime
 *   11. Dubins ground risk
 *   12. Dubins airspace risk
 *   13. Dubins flight-path-angle change
 *   14. Dubins runtime
 *   15. Dubins risk-evaluation runtime
 *   16. Selected landing-site identifier
 *   17. Highest-ranked landing-site identifier
 *
 * Returns:
 *   0 on success
 *  -1 on failure
 */
static int appendLogEntry(
    const char *logPath,
    const SearchProblem *problem,
    PathMetrics searchMetrics,
    PathMetrics dubinsMetrics,
    int selectedSiteIdent,
    int highestRankedSiteIdent) {

    FILE *logFile = fopen(logPath, "a");

    if (logFile == NULL) {
        perror("Unable to open result log");
        return -1;
    }

    const int writeStatus = fprintf(
        logFile,
        "%d,%.6f,%.6f,%.1f,%.1f,"
        "%.6f,%.6f,%.3f,%.0f,%.0f,"
        "%.6f,%.6f,%.3f,%.0f,%.0f,"
        "%d,%d\n",
        problem->exitFlag,
        problem->Initial.lat,
        problem->Initial.lon,
        problem->Initial.alt,
        problem->Initial.hdg,
        searchMetrics.groundRisk,
        searchMetrics.airspaceRisk,
        searchMetrics.dGamma,
        searchMetrics.totalRuntime,
        searchMetrics.riskRuntime,
        dubinsMetrics.groundRisk,
        dubinsMetrics.airspaceRisk,
        dubinsMetrics.dGamma,
        dubinsMetrics.totalRuntime,
        dubinsMetrics.riskRuntime,
        selectedSiteIdent,
        highestRankedSiteIdent);

    if (writeStatus < 0) {
        fprintf(
            stderr,
            "Unable to write result to %s\n",
            logPath);

        fclose(logFile);
        return -1;
    }

    if (fclose(logFile) != 0) {
        perror("Unable to close result log");
        return -1;
    }

    return 0;
}


/*
 * Return the identifier of the landing site selected by the planner.
 *
 * A return value of -1 indicates that no valid selected site is
 * available.
 */
static int getSelectedSiteIdent(const SearchProblem *problem) {
    if (problem == NULL ||
        problem->rankedSites == NULL ||
        problem->siteIndex < 0 ||
        problem->rankedSites[problem->siteIndex] == NULL) {
        return -1;
    }

    return problem->rankedSites[problem->siteIndex]->ident;
}


/*
 * Return the identifier of the highest-ranked landing site.
 *
 * A return value of -1 indicates that no ranked landing site is
 * available.
 */
static int getHighestRankedSiteIdent(const SearchProblem *problem) {
    if (problem == NULL ||
        problem->rankedSites == NULL ||
        problem->rankedSites[0] == NULL) {
        return -1;
    }

    return problem->rankedSites[0]->ident;
}


/*
 * Release all resources owned by this executable.
 *
 * Pointer-to-pointer arguments allow the function to clear each
 * caller-side pointer after releasing the corresponding resource.
 *
 * The alias check prevents a double free when minRiskDubins() returns
 * the same object that was supplied as the workspace.
 *
 */
static void cleanupResources(
    SearchProblem **problem,
    struct DubinsPath **workspace,
    struct DubinsPath **bestPath) {
    /*
     * Release the selected Dubins path only when it is different from
     * the workspace object.
     */
    if (bestPath != NULL && *bestPath != NULL) {
        if (workspace == NULL || *bestPath != *workspace) {
            free(*bestPath);
        }

        *bestPath = NULL;
    }

    /*
     * Release the Dubins workspace.
     */
    if (workspace != NULL && *workspace != NULL) {
        free(*workspace);
        *workspace = NULL;
    }

    /*
     * Release the planning problem.
     */
    if (problem != NULL && *problem != NULL) {
        freeProblem(*problem);
        *problem = NULL;
    }
}


int main(void) {
    /*
     * Mutable character arrays are used because the ACLM functions may
     * currently accept char pointers rather than const char pointers.
     *
     * These arrays are stack allocated and must not be freed.
     */
    char configPath[] = "aclm.cfg";
    char outputDirectory[] = "out/";
    char logPath[] = "out/log.csv";

    int programStatus = EXIT_SUCCESS;

    SearchProblem *problem = NULL;
    struct DubinsPath *dubinsWorkspace = NULL;
    struct DubinsPath *bestDubins = NULL;

    PathMetrics searchMetrics = unavailableMetrics();
    PathMetrics dubinsMetrics = unavailableMetrics();

    /*
     * Create the directory used for path files and the summary log.
     */
    if (ensureDirectoryExists(outputDirectory) != 0) {
        return EXIT_FAILURE;
    }

    /*
     * Allocate and initialize the planning problem.
     */
    problem = calloc(1, sizeof(*problem));

    if (problem == NULL) {
        perror("Unable to allocate SearchProblem");
        return EXIT_FAILURE;
    }

    /*
     * Read aclm.cfg, initialize the emergency scenario, rank candidate
     * landing sites, and execute the search-based contingency planner.
     */
    runEmergencyPlanning(
        problem,
        outputDirectory,
        configPath);

    printf("Exit flag: %d\n", problem->exitFlag);

    /*
     * Exit flags:
     *
     *  -1: No landing site is reachable.
     *
     *  -2: The emergency is initialized in prohibited airspace, or an
     *      invalid airspace-risk value was produced.
     */
    if (problem->exitFlag == -1 ||
        problem->exitFlag == -2) {

        if (appendLogEntry(
                logPath,
                problem,
                searchMetrics,
                dubinsMetrics,
                -1,
                -1) != 0) {
            programStatus = EXIT_FAILURE;
        }

        cleanupResources(
            &problem,
            &dubinsWorkspace,
            &bestDubins);

        return programStatus;
    }

    /*
     * Exit flags:
     *
     *   2: Search reached its state-expansion limit.
     *
     *  -3: Search exhausted the available landing-site candidates.
     *
     * In either case, attempt to recover a feasible solution using the
     * minimum-risk Dubins planner.
     */
    if (problem->exitFlag == 2 ||
        problem->exitFlag == -3) {

        double totalRuntime = 0.0;
        double riskRuntime = 0.0;

        dubinsWorkspace = calloc(
            1,
            sizeof(*dubinsWorkspace));

        if (dubinsWorkspace == NULL) {
            perror("Unable to allocate DubinsPath");

            cleanupResources(
                &problem,
                &dubinsWorkspace,
                &bestDubins);

            return EXIT_FAILURE;
        }

        bestDubins = minRiskDubins(
            dubinsWorkspace,
            problem,
            &totalRuntime,
            &riskRuntime);

        /*
         * No Dubins path was found, so the test case remains
         * unreachable.
         */
        if (bestDubins == NULL ||
            bestDubins->size == 0) {

            problem->exitFlag = -1;
        } else {
            /*
             * Evaluate the airspace risk of the fallback path.
             */
            airspaceRisk(
                problem,
                bestDubins,
                outputDirectory,
                1);

            /*
             * Ground-population risk is disabled for this test case.
             */
            if (!problem->w_gp) {
                bestDubins->gp = -1.0;
            }

            /*
             * A finite airspace-risk value indicates a valid fallback
             * path.
             */
            if (isfinite(bestDubins->ga)) {
                problem->exitFlag = 2;
            } else {
                problem->exitFlag = -2;
            }

            dubinsMetrics = (PathMetrics) {
                .groundRisk   = bestDubins->gp,
                .airspaceRisk = bestDubins->ga,
                .dGamma       = bestDubins->dgamma,
                .totalRuntime = totalRuntime,
                .riskRuntime  = riskRuntime
            };
        }

        if (appendLogEntry(
                logPath,
                problem,
                searchMetrics,
                dubinsMetrics,
                getSelectedSiteIdent(problem),
                getHighestRankedSiteIdent(problem)) != 0) {
            programStatus = EXIT_FAILURE;
        }

        /*
         * Write the fallback path only when a nonempty path exists.
         */
        if (bestDubins != NULL &&
            bestDubins->size > 0) {

            writeResults(
                problem,
                bestDubins,
                outputDirectory);
        }

        cleanupResources(
            &problem,
            &dubinsWorkspace,
            &bestDubins);

        return programStatus;
    }

    /*
     * The search-based planner found a candidate path. Evaluate its
     * ground risk if ground-population risk is enabled.
     */
    if (problem->w_gp) {
        problem->overflownPop = groundRisk(
            problem,
            NULL,
            outputDirectory,
            0);
    } else {
        problem->overflownPop = -1.0;
    }

    /*
     * Evaluate the airspace risk of the search-based path.
     */
    problem->airspaceOccup = airspaceRisk(
        problem,
        NULL,
        outputDirectory,
        0);

    searchMetrics = (PathMetrics) {
        .groundRisk   = problem->overflownPop,
        .airspaceRisk = problem->airspaceOccup,
        .dGamma       = problem->dGamma,
        .totalRuntime = problem->totalSearchRuntime,
        .riskRuntime  = problem->totalRiskRuntime
    };

    /*
     * Reject a search result whose airspace-risk calculation produced
     * an invalid value.
     */
    if (!isfinite(problem->airspaceOccup)) {
        problem->exitFlag = -2;

        fprintf(
            stderr,
            "Search path produced an invalid airspace-risk value.\n");

        if (appendLogEntry(
                logPath,
                problem,
                searchMetrics,
                dubinsMetrics,
                getSelectedSiteIdent(problem),
                getHighestRankedSiteIdent(problem)) != 0) {
            programStatus = EXIT_FAILURE;
        }

        cleanupResources(
            &problem,
            &dubinsWorkspace,
            &bestDubins);

        return programStatus;
    }

    /*
     * Compute a minimum-risk Dubins path for comparison with the
     * search-based solution.
     *
     * This path is used only for benchmarking.
     */
    dubinsWorkspace = calloc(
        1,
        sizeof(*dubinsWorkspace));

    if (dubinsWorkspace == NULL) {
        perror("Unable to allocate DubinsPath");

        cleanupResources(
            &problem,
            &dubinsWorkspace,
            &bestDubins);

        return EXIT_FAILURE;
    }

    double dubinsTotalRuntime = 0.0;
    double dubinsRiskRuntime = 0.0;

    bestDubins = minRiskDubins(
        dubinsWorkspace,
        problem,
        &dubinsTotalRuntime,
        &dubinsRiskRuntime);

    /*
     * The search-based path remains valid even if the Dubins benchmark
     * cannot produce a path.
     */
    if (bestDubins == NULL ||
        bestDubins->size == 0) {

        problem->exitFlag = 4;
    } else {
        /*
         * Evaluate the airspace risk of the Dubins benchmark path.
         */
        airspaceRisk(
            problem,
            bestDubins,
            outputDirectory,
            1);

        /*
         * Ground-population risk is disabled for this test case.
         */
        if (!problem->w_gp) {
            bestDubins->gp = -1.0;
        }

        /*
         * Mark the benchmark as invalid if its total risk is not finite.
         */
        if (!isfinite(bestDubins->risk)) {
            problem->exitFlag = 4;
        }

        dubinsMetrics = (PathMetrics) {
            .groundRisk   = bestDubins->gp,
            .airspaceRisk = bestDubins->ga,
            .dGamma       = bestDubins->dgamma,
            .totalRuntime = dubinsTotalRuntime,
            .riskRuntime  = dubinsRiskRuntime
        };
    }

    /*
     * Record the search and Dubins results in one CSV row.
     */
    if (appendLogEntry(
            logPath,
            problem,
            searchMetrics,
            dubinsMetrics,
            getSelectedSiteIdent(problem),
            getHighestRankedSiteIdent(problem)) != 0) {
        programStatus = EXIT_FAILURE;
    }

    /*
     * Write path-coordinate files when a valid Dubins benchmark path
     * is available.
     */
    if (bestDubins != NULL &&
        bestDubins->size > 0) {

        writeResults(
            problem,
            bestDubins,
            outputDirectory);
    }

    /*
     * Release all dynamically allocated resources before exiting.
     */
    cleanupResources(
        &problem,
        &dubinsWorkspace,
        &bestDubins);

    return programStatus;
}