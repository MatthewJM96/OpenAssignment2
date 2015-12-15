/**
* DeliveryOptimiser.cpp
* Takes a collection of delivery locations and optimises the placement of a hub to service those locations,
* as well as constructing a plan for delivery. If multiple hubs are desired, locations are separated into
* multiple regions.
*
* Author: Matthew Marshall
* Date Created: 03/12/2015
*/

#include <iostream>
#include <conio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <fstream>
#include <sstream>
#include <random>
#include <climits>
#include <stdarg.h>

/************************************************************************/
/* Helpers                                                              */
/************************************************************************/

#define PI 3.14159
#define EARTH_RADIUS 3958.75

/// Little helper to exit the program after a key press.
inline void exitProgram(size_t code) {
    std::cout << std::endl << "Press any key to exit..." << std::endl;
    _getch();
    exit(code);
}

/// Collection of helper function to act on strings.
namespace StringHelper {
    /// Splits a string by the provided delimiter.
    inline unsigned int split(const std::string& txt, std::vector<std::string>& strs, char ch)
    {
        unsigned int pos = txt.find(ch);
        unsigned int initialPos = 0;
        strs.clear();

        // Decompose statement
        while (pos != std::string::npos) {
            strs.push_back(txt.substr(initialPos, pos - initialPos));
            initialPos = pos + 1;

            pos = txt.find(ch, initialPos);
        }

        // Add the last one
        strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

        return strs.size();
    }

    /// Trims a string from the left.
    inline std::string& ltrim(std::string& s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    /// Trims a string from the right.
    inline std::string& rtrim(std::string& s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    /// Trims a string from both ends.
    inline std::string& trim(std::string& s) {
        return ltrim(rtrim(s));
    }
}

/// Provides the number of digits for an integer value.
template <typename T>
inline int iNumberOfDigits(T number)
{
    // Allows counting of integer digits of a floating-point or decimal type by static casting to largest signed integer type.
    intmax_t temp = static_cast<intmax_t>(std::round(number));
    int digits = 0;
    if (temp == 0) return 1;
    if (temp < 0) digits = 1;
    while (temp) {
        temp /= 10;
        digits++;
    }
    return digits;
}

/// Generates random values, either between min and max continuously or in steps.
inline double generateRandomNumber(double min, double max) {
    return min + (max * rand() / RAND_MAX);
}
inline double generateRandomNumber(double min, double max, unsigned int stepCount) {
    return min + (rand() % (stepCount + 1) * (1.0 / stepCount) * (max - min));
}

/// Grabs an integer from the user that is valid and in the range specified.
inline int getIntegerFromUser(int min, int max) {
    std::string line;
    int integer;
    while (std::getline(std::cin, line)) {
        std::stringstream lineStream(line);
        if (lineStream >> integer) {
            if (lineStream.eof()) {
                if (integer > max || integer < min) {
                    std::cout << "Please input an integer between " << min << " and " << max << "!" << std::endl;
                    continue;
                }
                break;
            }
        }
        std::cout << "Please input an integer!" << std::endl;
    }
    return integer;
}

/// Calculates the great circle distance between two points.
inline double greatCircleDistanceRad(double lat1, double lat2, double long1, double long2, double cosLat1) {
    double dLat = lat2 - lat1;
    double dLong = long2 - long1;

    double a = std::pow(std::sin(dLat / 2.0), 2.0) + cosLat1 * std::cos(lat2) * std::pow(std::sin(dLong / 2.0), 2.0);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));

    return c * EARTH_RADIUS;
}
inline double greatCircleDistanceRad(double lat1, double lat2, double long1, double long2) {
    double costLat1 = std::cos(lat1);

    return greatCircleDistanceRad(lat1, lat2, long1, long2, costLat1);
}
inline double greatCircleDistanceDeg(double lat1, double lat2, double long1, double long2) {
    double lat1Rad = lat1 / 360.0 * 2 * PI;
    double long1Rad = long1 / 360.0 * 2 * PI;
    double lat2Rad = lat2 / 360.0 * 2 * PI;
    double long2Rad = long2 / 360.0 * 2 * PI;

    return greatCircleDistanceRad(lat1Rad, lat2Rad, long1Rad, long2Rad);
}

/************************************************************************/
/* Model                                                                */
/************************************************************************/

/// Enumeration of location types.
enum class LocationType {
    CITY,
    TOWN
};

/// Information for a given location.
struct Location {
    std::string name;
    LocationType type;
    unsigned int population;
    double latitude;
    double longitude;
};

class LocationDataModel {
public:
    LocationDataModel() {}
    /// Ensure file is closed and locations cleared.
    ~LocationDataModel() {
        m_file.close();
    }

    /// Reads the specified file.
    bool readFile(char separator) {
        if (!m_file.is_open() && !openFile('i')) {
            return false;
        }
        clearLocations();
        readData(separator);
        return true;
    }

    // Setters
    /// Set the file path to.
    void setFilePath(const std::string& filePath) { m_filePath = filePath; }
    /// Set the location data.
    void setLocationData(const std::vector<Location>& locations) { m_locations = locations; }

    // Getters
    /// Get the currently set file path for data to be read from/written to.
    std::string getFilePath() const { return m_filePath; }
    /// Get the read location data.
    std::vector<Location> getLocationData() const { return m_locations; }
private:
    /// Opens the file at the currently set file path.
    /// Checks if the file is open and returns true on success, false on failure.
    bool openFile(char io = 'i') {
        // Open in correct mode.
        if (io == 'i') {
            m_file.open(m_filePath, std::ios::in);
        } else if (io == 'o') {
            m_file.open(m_filePath, std::ios::out);
        } else {
            return false;
        }
        if (m_file.is_open()) {
            return true;
        }
        return false;
    }

    /// Reads data from the currently open file.
    void readData(char separator) {
        // Read file lines into data points.
        std::string line;
        while (std::getline(m_file, line)) {
            // Skip if a comment.
            if (line.find('%') == 0) continue;
            
            Location row;
            std::vector<std::string> rowTemp;
            // Split a line by whitespace.
            StringHelper::split(line, rowTemp, separator);
            // Ensure row is of expected length.
            if (rowTemp.size() != 5) {
                std::cout << "Invalid file formatting. Could not read file." << std::endl;
                exitProgram(5);
            }
            // Trim the row elements of remaining whitespace.
            trimRow(rowTemp);
            // Populate the row with valid doubles.
            populateRow(rowTemp, row);
            // Add row to locations.
            m_locations.push_back(row);
        }
    }

    /// Trims an entire vector of strings of whitespace.
    std::vector<std::string>& trimRow(std::vector<std::string>& row) {
        for (auto& it : row) {
            StringHelper::trim(it);
        }
        return row;
    }

    /// Populates Location vector row with the contents of a given string vector row.
    void populateRow(const std::vector<std::string>& stringRow, Location& destRow) {
        destRow.name = stringRow[0];

        if (stringRow[1] == "City") {
            destRow.type = LocationType::CITY;
        } else if (stringRow[1] == "Town") {
            destRow.type = LocationType::TOWN;
        } else {
            std::cout << "Invalid value for one of the 'type' values." << std::endl;
            exitProgram(6);
        }

        try {
            destRow.population = stoi(stringRow[2]);
        } catch (std::invalid_argument) {
            std::cout << "Invalid value for one of the 'population' values." << std::endl;
            exitProgram(7);
        }

        try {
            destRow.latitude = stod(stringRow[3]);
        } catch (std::invalid_argument) {
            std::cout << "Invalid value for one of the 'latitude' values." << std::endl;
            exitProgram(7);
        }

        try {
            destRow.longitude = stod(stringRow[4]);
        } catch (std::invalid_argument) {
            std::cout << "Invalid value for one of the 'longitude' values." << std::endl;
            exitProgram(8);
        }
    }

    /// Clears the location vector.
    void clearLocations() {
        m_locations.clear();
    }

    /// File object.
    std::fstream m_file;
    /// File path string.
    std::string m_filePath = "GBplaces.csv";
    /// Collection of locations.
    std::vector<Location> m_locations;
};

/************************************************************************/
/* Hill Climb Algorithm                                                 */
/************************************************************************/

/// Hill climbing point structure.
struct HC_Point2D {
    double x;
    double y;
    double fitness;
};

template <typename... Args> ///< Allow fitness function pointers to have some statefulness without giving them real state.
class HillClimb2D {
public:
    HillClimb2D() {}
    /// Ensure function pointer is nullified.
    ~HillClimb2D() {
        m_fitness = nullptr;
    }

    /// Initialise a HillClimb2D object with given parameters.
    void init(double xMin, double yMin, double xMax, double yMax, double xStep, double yStep, double(*fitness)(double, double, Args...)) {
        m_xMin = xMin;  m_xMax = xMax;  m_xStep = xStep;
        m_yMin = yMin;  m_yMax = yMax;  m_yStep = yStep;
        m_fitness = fitness;
    }

    /// Run hill climb algorithm runCount times, choosing best fitness.
    HC_Point2D run(unsigned int runCount, Args... args) {
        HC_Point2D bestMaxima = { 0.0, 0.0, std::numeric_limits<double>::lowest() };
        for (size_t i = 0; i < runCount; ++i) {
            HC_Point2D currentMaxima = climbHill(args...);
            if (currentMaxima.fitness > bestMaxima.fitness) {
                bestMaxima = currentMaxima;
            }
        }
        return bestMaxima;
    }
private:
    /// Run hill climb algorithm returning the achieved maxima's point.
    HC_Point2D climbHill(Args... args) {
        // Generate a random starting point.
        double x = generateRandomNumber(m_xMin, m_xMax, static_cast<unsigned int>((m_xMax - m_xMin) / m_xStep));
        double y = generateRandomNumber(m_yMin, m_yMax, static_cast<unsigned int>((m_yMax - m_yMin) / m_yStep));

        // Get starting point's fitness.
        double currentFitness = m_fitness(x, y, args...);
        double oldFitness;

        do {
            oldFitness = currentFitness;

            // Loop over adjacent (diagonally, vertically and horizontally) cells, checking for better fitnesses than the current.
            int dx = 0;
            int dy = 0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    // Skip self.
                    if (i == 0 && j == 0) continue;

                    double xCurr = x + m_xStep * i;
                    double yCurr = y + m_yStep * j;

                    // Skip out of bounds cells.
                    if (xCurr > m_xMax || xCurr < m_xMin ||
                        yCurr > m_yMax || yCurr < m_yMin) continue;

                    // If new fitness is better than current, store the direction of the better cell, and set current fitness to the new value.
                    double newFitness = m_fitness(xCurr, yCurr, args...);
                    if (newFitness >= currentFitness) {
                        dx = i;
                        dy = j;
                        currentFitness = newFitness;
                    }
                }
            }

            // Update "pivot" cell to be the best found.
            x += m_xStep * dx;
            y += m_yStep * dy;
        } while (currentFitness > oldFitness);

        return HC_Point2D{ x, y, oldFitness };
    }

    double m_xMin, m_xMax, m_xStep, m_yMin, m_yMax, m_yStep;
    double(*m_fitness) (double, double, Args...);
};

/************************************************************************/
/* K-Means Clustering Algorithm                                         */
/************************************************************************/

struct KMC_Point2D {
    double x, y, f;
    unsigned int group;
    unsigned int foreignId;
};

template <typename... Args>
class KMeansClustering2D {
public:
    KMeansClustering2D() {}
    /// Ensure function pointer is nullified.
    ~KMeansClustering2D() {
        m_distanceToClusterCenter = nullptr;
    }

    /// Initialise a Lloyds2D object with given parameters.
    void init(const std::vector<KMC_Point2D>& points, unsigned int clusterCount, double(*distanceToClusterCenter) (KMC_Point2D, KMC_Point2D, Args...)) {
        m_points = points;
        m_clusterCount = clusterCount;
        m_distanceToClusterCenter = distanceToClusterCenter;
    }

    // TODO(Matthew): Implement a multiple run count option to reduce changes of local maxima not == global maxima.
    /// Run Lloyd's algorithm, returning clusters.
    std::vector<std::vector<KMC_Point2D>> run(size_t accecptableChangeRate, size_t maxIterations, bool& maxIterationCountReached, Args... args) {
        if (m_clusterCount == 0) throw;
        if (m_clusterCount == 1) {
            std::vector<std::vector<KMC_Point2D>> clusters;
            clusters.push_back(m_points);
            return clusters;
        }
        return kMeansClustering(accecptableChangeRate, maxIterations, maxIterationCountReached, args...);
    }
private:
    /// Performs Lloyd's algorithm, returning the collection of clusters constructed.
    std::vector<std::vector<KMC_Point2D>> kMeansClustering(size_t accecptableChangeRate, size_t maxIterations, bool& maxIterationCountReached, Args... args) {
        std::vector<KMC_Point2D> clusterCenters;
        clusterCenters.resize(m_clusterCount);

        // Prepare cluster centers using the k-means++ algorithm.
        kpp(clusterCenters, args...);

        size_t iterations = 1;
        unsigned int changed;
        do {
            changed = 0;

            std::vector<KMC_Point2D> newClusterCenters;
            std::vector<int> newClusterPointCounts;
            newClusterCenters.resize(m_clusterCount, {});
            newClusterPointCounts.resize(m_clusterCount, 0);

            for (size_t i = 0; i < m_points.size(); ++i) {
                // Associate point with nearest cluster center.
                int closestClusterCenter = nearestClusterCenter(m_points[i], clusterCenters);
                if (closestClusterCenter != m_points[i].group) {
                    ++changed;
                    m_points[i].group = closestClusterCenter;
                }
                // Add position to new cluster position for calculating centroid of associated cluster points.
                newClusterCenters[m_points[i].group].x += m_points[i].x;
                newClusterCenters[m_points[i].group].y += m_points[i].y;
                newClusterCenters[m_points[i].group].f += m_points[i].f;
                ++newClusterPointCounts[m_points[i].group];
            }

            // Calculate centroids and set them as new cluster centers.
            for (size_t i = 0; i < m_clusterCount; ++i) {
                newClusterCenters[i].x /= newClusterPointCounts[i];
                newClusterCenters[i].y /= newClusterPointCounts[i];
            }
            clusterCenters = newClusterCenters;
            if (++iterations > maxIterations) {
                maxIterationCountReached = true;
                break;
            }
        } while (changed > accecptableChangeRate);

        std::vector<std::vector<KMC_Point2D>> clusters;
        clusters.resize(m_clusterCount);
        for (auto& point : m_points) {
            clusters[point.group].push_back(point);
        }

        return clusters;
    }

    /// k-means++ algorithm for specifying beginning cluster centers (results in faster and better clusters).
    void kpp(std::vector<KMC_Point2D>& clusterCenters, Args... args) {
        size_t preparedClusterCentersCount = 0;

        // Randomly select the first cluster center as an existing point.
        clusterCenters[0] = m_points[rand() % m_points.size()];
        clusterCenters[0].f = 0;
        ++preparedClusterCentersCount;

        // Repeat contents for each cluster center.
        for (size_t i = 1; i < m_clusterCount; ++i) {
            std::vector<double> distancesSquared;
            // Calculate distances to nearest of prepared clusters.
            for (auto& point : m_points) {
                double shortestDistance = std::numeric_limits<double>::max();
                for (size_t j = 0; j < preparedClusterCentersCount; ++j) {
                    double dist = m_distanceToClusterCenter(point, clusterCenters[j], args...);
                    if (dist < shortestDistance) {
                        shortestDistance = dist;
                    }
                }
                distancesSquared.push_back(std::pow(shortestDistance, 2.0));
            }

            // Calculate total distances.
            double totalDistancesSquared = 0.0;
            for (auto& dist : distancesSquared) {
                totalDistancesSquared += dist;
            }

            // Randomly select a point as a new cluster center using min existing cluster center distance as weighting.
            double randVal = generateRandomNumber(0.0, 1.0);
            for (size_t j = 0; j < distancesSquared.size(); ++j) {
                double distNormalised = distancesSquared[j] / totalDistancesSquared;
                if (randVal < distNormalised) {
                    clusterCenters[i] = m_points[j];
                    clusterCenters[i].f = 0;
                    ++preparedClusterCentersCount;
                    break;
                }
                randVal -= distNormalised;
            }
        }
    }

    // Two separate functions for index and distance of nearest cluster center for performance-sake.
    /// Returns index of nearest cluster center to given point.
    int nearestClusterCenter(KMC_Point2D point, std::vector<KMC_Point2D> clusterCenters, Args... args) {
        int index = -1;
        double shortestDistance = std::numeric_limits<double>::max();
        for (size_t i = 0; i < clusterCenters.size(); ++i) {
            double dist = m_distanceToClusterCenter(point, clusterCenters[i], args...);
            if (dist < shortestDistance) {
                shortestDistance = dist;
                index = i;
            }
        }
        return index;
    }
    /// Returns distance of nearest cluster center to given point.
    double distanceToNearestClusterCenter(KMC_Point2D point, std::vector<KMC_Point2D> clusterCenters, Args... args) {
        double shortestDistance = std::numeric_limits<double>::max();
        for (size_t i = 0; i < clusterCenters.size(); ++i) {
            double dist = m_distanceToClusterCenter(point, clusterCenters[i], args...);
            if (dist < shortestDistance) {
                shortestDistance = dist;
            }
        }
        return shortestDistance;
    }

    unsigned int m_clusterCount;
    std::vector<KMC_Point2D> m_points;
    double(*m_distanceToClusterCenter) (KMC_Point2D, KMC_Point2D, Args...);
};

/************************************************************************/
/* View                                                                 */
/************************************************************************/

class GenericTextScreen {
public:
    /// Render content to screen.
    virtual void render() {
        printBuffer();
    }

    /// Flushes buffer to screen.
    virtual void flush() {
        printBuffer();
        m_buffer.str(std::string());
    }

    /// Adds a string to the buffer for the screen.
    virtual GenericTextScreen* addToBuffer(const std::string& message) {
        m_buffer << message;
        return this;
    }
protected:
    /// Prints buffer to screen.
    virtual void printBuffer() const {
        std::cout << m_buffer.str();
    }

    std::stringstream m_buffer;
};

class HubLocatorScreen : public GenericTextScreen {
public:
    /// Include generic addToBuffer functions.
    using GenericTextScreen::addToBuffer;
    /// Add a location to buffer for the screen.
    virtual HubLocatorScreen* addToBuffer(const Location& location) {
        m_buffer << "    Name        -  " << location.name << "\n";
        m_buffer << "    Type        -  " << (location.type == LocationType::CITY ? "City" : "Town") << "\n";
        m_buffer << "    Population  -  " << location.population << "\n";
        m_buffer << "    Latitude    -  " << location.latitude << "\n";
        m_buffer << "    Longitude   -  " << location.longitude << "\n";
        return this;
    }
    /// Add a location to buffer for the screen.
    virtual HubLocatorScreen* addToBuffer(const std::vector<Location>& region) {
        m_buffer << "    Locations:\n";
        unsigned int population = 0;
        for (auto& location : region) {
            population += location.population;
            m_buffer << "        - " << location.name << "\n";
        }
        m_buffer << "    Population  - " << population << "\n";
        return this;
    }
    /// Add a hub point to buffer for the screen.
    virtual HubLocatorScreen* addToBuffer(const HC_Point2D& hub) {
        m_buffer << "    Latitude    -  " << hub.x << "\n";
        m_buffer << "    Longitude   -  " << hub.y << "\n";
        return this;
    }
};

/************************************************************************/
/* Controller                                                           */
/************************************************************************/

class HubPlacer {
public:
    HubPlacer() {}
    ~HubPlacer() {}

    void run(unsigned int hubCount, std::string filePath) {
        LocationDataModel model;

        m_view.addToBuffer("Reading file...\n")->flush();

        // Read in file.
        model.setFilePath(filePath);
        model.readFile(',');
        std::vector<Location> locations = model.getLocationData();

        m_view.addToBuffer("File read.\n\nPreparing region constructor...\n")->flush();

        // Construct vector of Lloyd's algorithm points.
        std::vector<KMC_Point2D> lPoints;
        for (size_t i = 0; i < locations.size(); ++i) {
            KMC_Point2D lPoint = {
                locations[i].latitude,
                locations[i].longitude,
                static_cast<double>(locations[i].population),
                0,
                i,
            };
            lPoints.push_back(lPoint);
        }

        m_view.addToBuffer("Constructing regions...")->flush();

        // Run Lloyd's algorithm to obtain clusters of locations.
        KMeansClustering2D<> lManager;
        lManager.init(lPoints, hubCount,
            [](KMC_Point2D point, KMC_Point2D center) -> double {
                return greatCircleDistanceDeg(point.x, center.x, point.y, center.y) * std::sqrt(center.f + point.f);
            }
        );
        bool maxIterationCountReached = false;
        std::vector<std::vector<KMC_Point2D>> clusters = lManager.run(0, 10000, maxIterationCountReached);

        // Convert L_Point2D clusters back to collections of locations into regions.
        std::vector<std::vector<Location>> regions;
        for (auto& cluster : clusters) {
            if (cluster.empty()) {
                continue;
            }
            std::vector<Location> region;
            for (auto& point : cluster) {
                region.push_back(locations[point.foreignId]);
            }
            regions.push_back(region);
            m_view.addToBuffer("\nRegion #" + std::to_string(regions.size()) + " constructed.\n");
            m_view.addToBuffer(regions.back())->flush();
        }

        m_view.addToBuffer("\n\nRegions constructed.\n\nPreparing hub locators.\n\n")->flush();

        // Get the optimised hub locations for each region.
        std::vector<HC_Point2D> hubLocations;
        for (auto& region : regions) {
            HillClimb2D<std::vector<Location>> hcManager;
            hcManager.init(50.0, -5.0, 58.0, 2.0, 0.005, 0.005,
                [](double lat, double lon, std::vector<Location> locations) -> double {
                    double fitness = 0.0;

                    // Minimise operations.
                    double selfLatRad = lat / 360.0 * 2 * PI;
                    double selfLongRad = lon / 360.0 * 2 * PI;
                    double cosSelfLatRad = std::cos(selfLatRad); 

                    // Calculate distance to each location from the specified point and weight with population of a given location.
                    for (auto& it : locations) {
                        double locLatRad = it.latitude / 360.0 * 2 * PI;
                        double locLongRad = it.longitude / 360.0 * 2 * PI;

                        fitness -= greatCircleDistanceRad(selfLatRad, locLatRad, selfLongRad, locLongRad, cosSelfLatRad) * it.population;
                    }

                    return fitness;
                }
            );

            m_view.addToBuffer("Hub locator prepared.\nCalculating hub location...\n")->flush();

            hubLocations.push_back(hcManager.run(1, region));

            m_view.addToBuffer("Hub location calculated.\n\n");
        }

        if (maxIterationCountReached) {
            m_view.addToBuffer("\nWARNING: Max iteration count reached in region calculation, displayed hub locations may be sub-optimal.\n");
        }

        if (hubLocations.size() != hubCount) {
            m_view.addToBuffer("\nWARNING: The number of regions constructed did not match the number requested.\n");
        }

        m_view.addToBuffer("\nHub locations calculated.\n\n");

        for (size_t i = 0; i < hubLocations.size(); ++i) {
            m_view.addToBuffer("Hub #" + std::to_string(i + 1) + ":\n");
            m_view.addToBuffer(hubLocations[i]);
            m_view.addToBuffer("\n");
        }

        m_view.flush();
    }
private:
    HubLocatorScreen m_view;
};

/************************************************************************/
/* Program                                                              */
/************************************************************************/

int main(int argc, char** argv) {
    // Seed random generator.
    srand(static_cast<unsigned int>(time(nullptr)));

    std::cout << "Note: In the event that regions constructed do not allow for the specified number of hubs, hubs will only be placed for regions constructed." << std::endl << std::endl;
    std::cout << "How many hubs would you like to place? (0-100):" << std::endl;
    int hubCount = getIntegerFromUser(0, 100);
    std::cout << std::endl;

    HubPlacer hp;
    hp.run(hubCount, "GBplaces.csv");

    std::cout << std::endl << "Press any key to exit..." << std::endl;
    _getch();

    return 0;
}