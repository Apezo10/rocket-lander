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
    static constexpr double length = 124.4; //meters

    //Mass will chnage once thrust is implemeted
    static constexpr double mass = 100000; // kg
    double frontArea = pi * diameter * diameter / 4.0;

    struct xyz {
        double x{};
        double y{};
        double z{};
    };

    struct inertia {
        double Ixx{};
        double Iyy{};
        double Izz{};
    };

    struct quaternion {
        double w{1.0};
        double x{0.0};
        double y{0.0};
        double z{0.0};
    };

    double dt = 0.001;
    double grav{};

    xyz initialPos;
    xyz currentPos;
    xyz initialVel;
    xyz currentVel;

    xyz angularVelocity{};
    xyz angularAcceleration{};
    xyz torque{};


    double calculateInertiaPY() {
        return (1.0 / 12.0) * mass * (length * length) + (1.0/4.0) * mass * ((diameter/2.0) * (diameter/2.0));
    }

    double calculateInertiaR() {
        return (1.0/2.0) * mass * ((diameter/2.0) * (diameter/2.0));
    }

    quaternion qMultiply(const quaternion& a, const quaternion& b) {
        return {
                a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,

                a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,

                a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,

                a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
            };
    }

    //initializes orientation in the upward position
    quaternion orientation{std::sqrt(0.5),0.0,-std::sqrt(0.5),0.0};

    //innitalize moment of inertia in roll, pitch, and yaw
    inertia rocketInertia{
        calculateInertiaR(),
        calculateInertiaPY(),
        calculateInertiaPY()
    };

    public:
    std::vector<xyz> pos{};
    std::vector<xyz> vel{};
    std::vector<xyz> angularVel{};
    std::vector<quaternion> orientations{};

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

        double dragFactor = -0.5 * rho * Cd * frontArea * speed;

        drag.x = dragFactor * currentVel.x;
        drag.y = dragFactor * currentVel.y;
        drag.z = dragFactor * currentVel.z;

        return drag;
    }

    // Moment of inertia for pitch and yaw axis
    double getInertiaPY() {
        return calculateInertiaPY();
    }

    //Moment of inertia for roll axis
    double getInertiaR() {
        return calculateInertiaR();
    }

    void calculateAngularAccel() {
        double p = angularVelocity.x;
        double q = angularVelocity.y;
        double r = angularVelocity.z;

        angularAcceleration.x = (torque.x + 
                                (rocketInertia.Iyy - rocketInertia.Izz)
                                 * q * r) / rocketInertia.Ixx;
                                
        angularAcceleration.y = (torque.y + 
                                (rocketInertia.Izz - rocketInertia.Ixx)
                                 * r * p) / rocketInertia.Iyy;
                                 
        angularAcceleration.z = (torque.z + 
                                (rocketInertia.Ixx - rocketInertia.Iyy)
                                 * p * q) / rocketInertia.Izz;
    }

    void updateAngularrVel(double dt) {
        angularVelocity.x += angularAcceleration.x * dt;
        angularVelocity.y += angularAcceleration.y * dt;
        angularVelocity.z += angularAcceleration.z * dt;
    }

    void updateOrientation(double dt) {

        quaternion omegaQuat{0.0, angularVelocity.x, angularVelocity.y, angularVelocity.z};

        quaternion qDot = qMultiply(orientation, omegaQuat);

        orientation.w += 0.5 * qDot.w * dt;
        orientation.x += 0.5 * qDot.x * dt;
        orientation.y += 0.5 * qDot.y * dt;
        orientation.z += 0.5 * qDot.z * dt;

        //normalize
        double norm = std::sqrt(
                orientation.w * orientation.w +
                orientation.x * orientation.x +
                orientation.y * orientation.y +
                orientation.z * orientation.z
            );

        orientation.w /= norm;
        orientation.x /= norm;
        orientation.y /= norm;
        orientation.z /= norm;
    }


    void updateState(const earth& e) {

        angularVel.push_back(angularVelocity);
        orientations.push_back(orientation);

        while (true) {

            //  1) Get gravity and air density as function of altitude
            double gravity = e.getGravity(currentPos.z);
            double rho = e.getAirDensity(currentPos.z);

            //  2) Get drag as function of density (USE 0.8 as Cd for now until we get accurate dynamic estimates)
            xyz drag = getDrag(rho, 0.8);

            // 3) Get accel
            xyz accel{};
            accel.x = drag.x / mass;
            accel.y = drag.y / mass;
            accel.z = drag.z / mass - gravity;
            calculateAngularAccel();

            //  4) Update vel
            currentVel.x += accel.x * dt;
            currentVel.y += accel.y * dt;
            currentVel.z += accel.z * dt;
            updateAngularrVel(dt);

            //  5) Update pos
            currentPos.x += currentVel.x * dt;
            currentPos.y += currentVel.y * dt;
            currentPos.z += currentVel.z * dt;
            updateOrientation(dt);

            //  6) store updated state
            vel.push_back(currentVel);
            pos.push_back(currentPos);
            angularVel.push_back(angularVelocity);
            orientations.push_back(orientation);

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

        file << "time,x,y,z,vx,vy,vz,p,q,r,qw,qx,qy,qz\n";

        //Iterates through pos and vel every second and writes to csv
        auto writeRow = [&](size_t i) {
            double time = i*dt;

            file << time << ","
             << pos[i].x << ","
             << pos[i].y << ","
             << pos[i].z << ","

             << vel[i].x << ","
             << vel[i].y << ","
             << vel[i].z << ","
             
             << angularVel[i].x << ","
             << angularVel[i].y << ","
             << angularVel[i].z << ","
             
             << orientations[i].w << ","
             << orientations[i].x << ","
             << orientations[i].y << ","
             << orientations[i].z << "\n";
        };

         // Export every second
        for (size_t i = 0; i < pos.size(); i += 1000) {
            writeRow(i);
        }

        // Export final state
        size_t last = pos.size() - 1;

        if (last % 1000 != 0) {
            writeRow(last);
        }

        file.close();

        std::cout << "CSV File exported successfully!\n";
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