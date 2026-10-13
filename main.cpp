#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>



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
    struct xyz {
        double x{};
        double y{};
        double z{};
    };

    double dx = 0.001;
    double grav{};

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

    void updateState(const earth& e) {

        //The for loop iterates through the timestep so this number determines the accuracy of the model
        for (size_t i = 0; i < 1000; i++) {

            double gravity = e.getGravity(currentPos.z);

            //Update position by integrating velocity in that direction
            currentPos.x += currentVel.x * dx;
            currentPos.y += currentVel.y * dx;
            currentPos.z += currentVel.z * dx - 0.5 * gravity * dx * dx;

            //Update velocity
            currentVel.z -= gravity * dx;
 
            //Store updated trajectory
            pos.push_back(currentPos);
            vel.push_back(currentVel);
        }
    }
};

int main() {
    rocket r;
    earth e;

    r.getInitialPos();
    r.getInitialVel();

    r.updateState(e);

    //current position
    for (int i = 0; i<r.pos.size(); i+=50) {
        std::cout << std::fixed << std::setprecision(3);

        std::cout <<  "X pos: " << r.pos[i].x
         << "\t" << "Y pos: " << r.pos[i].y 
         << "\t" << "Z pos: " << r.pos[i].z << "\n";
    }

    return 0;
}

//--------------------------
//NEXT STEP: Implemet aerodynamic drag and such other forces acting on the rocket. Also implement the surface of the planet so the rocket can actually hit the ground
//---------------------------------------------------------