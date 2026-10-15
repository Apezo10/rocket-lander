#include <iostream>
#include <cmath>
#include <vector>

class rocket {
    private:

    struct pos{
        double x{};
        double y{};
        double z{};
    };

    struct vel{
        double xVel{};
        double yVel{};
        double zVel{};
    };


    pos x;
    vel x_dot;

    std::vector<pos> rocketPos;
    std::vector<vel> rocketVel;

    public: 
};

