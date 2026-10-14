#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>



class earth {
    private:
    static constexpr double g = 9.80665;
    static constexpr double radius = 6371000.0;
    static constexpr double rho0 = 1.225;
    static constexpr double scaleHeight = 8500.0;

    public:
    double getGravity(double altitude) const {
        return g * std::pow(radius / (radius + altitude), 2);
    }

    double getAirDensity(double altitude) const {
         return rho0 * std::exp(-std::max(0.0, altitude) / scaleHeight);
    }
};


class rocket {
    private:
    static constexpr double diameter = 9.0;
    static constexpr double pi = 3.141592653589793;
    double area = pi * diameter * diameter / 4.0;

    struct xyz {
        double x{};
        double y{};
        double z{};
    };

    double dx = 0.001;
    double grav{};
    //Dry mass of startship booster in kg
    double dryMass = 180000;

    xyz initialPos;
    xyz currentPos;
    xyz initialVel;
    xyz currentVel;

    public:
    std::vector<xyz> pos{};
    std::vector<xyz> vel{};

    void getInitialPos() {
        std::cout << "What is the rockets x-coordniate (meters): ";
        std::cin >> initialPos.x;
        std::cout << "What is the rockets y-coordniate (meters): ";
        std::cin >> initialPos.y;
        std::cout << "What is the rockets z-coordniate (meters): ";
        std::cin >> initialPos.z;

        currentPos = initialPos;
        pos.push_back(initialPos);
    }

    void getInitialVel() {
        std::cout << "What is the rockets x-vel (m/s): ";
        std::cin >> initialVel.x;
        std::cout << "What is the rockets y-vel (m/s): ";
        std::cin >> initialVel.y;
        std::cout << "What is the rockets z-vel (m/s): ";
        std::cin >> initialVel.z;

        currentVel = initialVel;
        vel.push_back(initialVel);
    }


    //This calculates drag as a vector
    auto getDrag(double rho, double Cd) {
        xyz drag{};

        double speed = sqrt(currentVel.x * currentVel.x +
                            currentVel.y * currentVel.y +
                            currentVel.z * currentVel.z);

        double dragFactor = -0.5 * rho * Cd * area * speed;

        drag.x = dragFactor * currentVel.x;
        drag.y = dragFactor * currentVel.y;
        drag.z = dragFactor * currentVel.z;

        return drag;
    }

    void updateState(const earth& e) {

        //The for loop iterates through the timestep so this number determines the accuracy of the model
        while (true) {

            //  1) Get gravity and air density as function of altitude
            double gravity = e.getGravity(currentPos.z);
            double rho = e.getAirDensity(currentPos.z);

            //  2) Get drag as function of density (USE 0.8 as Cd for now until we get accurate dynamic estimates)
            xyz drag = getDrag(rho, 0.8);

            // 3) Get accel
            xyz accel{};
            accel.x = drag.x / dryMass;
            accel.y = drag.y / dryMass;
            accel.z = drag.z / dryMass - gravity;

            //  4) Update vel
            currentVel.x += accel.x * dx;
            currentVel.y += accel.y * dx;
            currentVel.z += accel.z * dx;

            //  5) Update pos
            currentPos.x += currentVel.x * dx;
            currentPos.y += currentVel.y * dx;
            currentPos.z += currentVel.z * dx;

            //  6) store updated state
            vel.push_back(currentVel);
            pos.push_back(currentPos);

            //Run sim until altitude is 0
            if (currentPos.z < 0.001) {
                currentPos.z = 0;

                pos.back() = currentPos;
                break;
            }
        }
    }

    void exportCSV(const std::string& filename) const {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error opening CSV file\n";
            return;
        }

        file << "time,x,y,z,vx,vy,vz\n";

        //Iterates through pos and vel every second and writes to csv
        for (size_t i = 0; i<pos.size(); i+=1000) {
            double time = i*dx;

            file << time << ","
             << pos[i].x << ","
             << pos[i].y << ","
             << pos[i].z << ","
             << vel[i].x << ","
             << vel[i].y << ","
             << vel[i].z << "\n";
        }

        // Always include the final state
        size_t last = pos.size() - 1;

        if (last % 1000 != 0) {
        file << last * dx << ","
             << pos[last].x << ","
             << pos[last].y << ","
             << pos[last].z << ","
             << vel[last].x << ","
             << vel[last].y << ","
             << vel[last].z << "\n";
        }

        file.close();

        std::cout << "CSV File exported successfully!";
    }
};

int main() {
    rocket r;
    earth e;

    r.getInitialPos();
    r.getInitialVel();

    r.updateState(e);

    //current position
    for (int i = 0; i<r.pos.size(); i+=100) {
        std::cout << std::fixed << std::setprecision(3);

        std::cout <<  "X pos: " << r.pos[i].x
         << "\t" << "Y pos: " << r.pos[i].y 
         << "\t" << "Z pos: " << r.pos[i].z 
         << "\t\t" << "X vel: " << r.vel[i].x 
         << "\t" << "Y vel: " << r.vel[i].y 
         << "\t" << "Z vel: " << r.vel[i].z << "\n";
    }

    r.exportCSV("Trajectory.csv");
    std::system("python handleCSV.py");

    return 0;
}

//---------------------------------------------------------
//NEXT STEP: Implement and api call to open mateo for accurate weather measurements based on lattitude and longitude up to 30 km. 
//           Use this to calculate accurate air desity for more precice drag measuremeents and wind measurements. For each api request,
//           ask for a vertical set of data based on varying airpressures to get approximate altitudes. Use this metho to request 
//           4 vertical regions which i can use to interpolate airpressure and wind based on rockets current xy (lat/long) pos. Only
//           request another call if the xy pos of the rocket moves beyond the grid formed by the 4 vertical regions. This way, we are
//           consistently using accurate data across the calculations.

// Does not assume a change in angle of attack for now (trjectory planning will come in the future)
//---------------------------------------------------------