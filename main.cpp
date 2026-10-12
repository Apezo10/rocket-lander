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

    xyz initialPos;
    xyz currentPos;
    xyz initialVel;
    xyz currentVel;

    public:
    std::vector<xyz> pos0{};
    std::vector<xyz> vel0{};

    void getInitialPos() {
        std::cout << "\nWhat is the rockets x-coordniate (meters): ";
        std::cin >> initialPos.x;
        std::cout << "\nWhat is the rockets y-coordniate (meters): ";
        std::cin >> initialPos.y;
        std::cout << "\nWhat is the rockets z-coordniate (meters): ";
        std::cin >> initialPos.z;

        pos0.push_back(initialPos);
    }

    void getInitialVel() {
        std::cout << "\nWhat is the rockets x-vel (m/s): ";
        std::cin >> initialVel.x;
        std::cout << "\nWhat is the rockets y-vel (m/s): ";
        std::cin >> initialVel.y;
        std::cout << "\nWhat is the rockets z-vel (m/s): ";
        std::cin >> initialVel.z;

        vel0.push_back(initialVel);
    }

    double getGravity(const std::vector<xyz>& positions) {
        float grav{};

        grav = g*((earthRadius/(earthRadius + positions[0].z))*(earthRadius/(earthRadius + positions[0].z)));
        return grav;
    }

};

int main() {
    rocket r;
    r.getInitialPos();

    std::cout << "grav: " << r.getGravity(r.pos0);

    return 0;
}
