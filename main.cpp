#include <iostream>
#include <cmath>
#include <vector>

static constexpr double g = 9.80665;
static constexpr double earthRadius = 6371000.0;



class rocket {
    private:
    struct xyz {
        float x{};
        float y{};
        float z{};
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

        pos.push_back(initialPos);
    }

    void getInitialVel() {
        std::cout << "What is the rockets x-vel (m/s): ";
        std::cin >> initialVel.x;
        std::cout << "What is the rockets y-vel (m/s): ";
        std::cin >> initialVel.y;
        std::cout << "What is the rockets z-vel (m/s): ";
        std::cin >> initialVel.z;

        vel.push_back(initialVel);
    }

    void updateState(const std::vector<xyz>& position, const std::vector<xyz>& velocity) {

        //The for loop iterates through the timestep so this number determines the accuracy of the model
        for (size_t i = 0; i < 1000; i++) {
            double gravity = g * std::pow(earthRadius / (earthRadius + pos[i].z), 2);

            currentPos.x = position[i].x + velocity[i].x*dx;
            currentPos.y = position[i].y + velocity[i].y*dx;
            currentPos.z = position[i].z + velocity[i].z*dx - 0.5 * gravity * dx * dx;

            currentVel.x = velocity[i].x;
            currentVel.y = velocity[i].y;
            currentVel.z = velocity[i].z - gravity * dx;
 
            pos.push_back(currentPos);
            vel.push_back(currentVel);
        }
    }
};

int main() {
    rocket r;
    r.getInitialPos();
    r.getInitialVel();
    r.updateState(r.pos, r.vel);


    //Prints the position of the rocket at 10 Hz
    for (size_t i=0; i<r.pos.size(); i+=100) {
        std::cout << "\nPosition is:" << "(" 
            << r.pos[i].x 
            << "," << r.pos[i].y
            << "," << r.pos[i].z << ")" 
            << "\n";    
        }


    return 0;
}
