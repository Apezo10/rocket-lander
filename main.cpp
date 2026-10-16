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

    double CoM = length/2.0;
    double CoP = 40.0; // fixed for now

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

    struct aero {
        xyz force{};
        xyz torque{};
    };

    //Implement a helper function for computing cross product
    xyz crossProduct(const xyz& a, const xyz& b) const {

        return {
            a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x
        };
    }

    double dt = 0.001;

    xyz initialPos;
    xyz currentPos;
    xyz initialVel;
    xyz currentVel;

    xyz angularVelocity{};
    xyz angularAcceleration{};
    xyz torque{};


    //Moment of inertia for pitch and yaw
    double calculateInertiaPY() {
        return (1.0 / 12.0) * mass * (length * length) + (1.0/4.0) * mass * ((diameter/2.0) * (diameter/2.0));
    }

    //Moment of inertia for roll
    double calculateInertiaR() {
        return (1.0/2.0) * mass * ((diameter/2.0) * (diameter/2.0));
    }

    //Implement a helper function for quaternion multiplication
    quaternion qMultiply(const quaternion& a, const quaternion& b) const {
        return {
                a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,

                a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,

                a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,

                a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
            };
    }

    //initializes orientation in the upward position
    quaternion orientation{
        std::sqrt(0.5), 0.0, -std::sqrt(0.5), 0.0
    };

    //innitalize moment of inertia in roll, pitch, and yaw
    inertia rocketInertia{
        calculateInertiaR(),
        calculateInertiaPY(),
        calculateInertiaPY()
    };

    xyz cpOffset {
        CoP - CoM,
        0.0,
        0.0
        };

    // Temporary aerodynamic normal-force coefficient slope
    static constexpr double normalForceSlope = 2.0;

    //Helper functions
    //--------------------------------------------------------------

    // 1) Rotate a vector using quaternions
    xyz rotateVector(const xyz& v, const quaternion& q) const {
        quaternion vectorQuat{0.0, v.x, v.y, v.z};

        quaternion qConjugate{
            q.w,
            -q.x,
            -q.y,
            -q.z
        };

        quaternion result = qMultiply(qMultiply(q, vectorQuat),qConjugate);

        return {
            result.x,
            result.y,
            result.z
        };
    };

    // 2) Convert a vector from body coordinates to world coordinates
    xyz bodyToWorld(const xyz& v) const {
        return rotateVector(v, orientation);
    };

    // 3) Convert a vector from world coordinates to body coordinates
    xyz worldToBody(const xyz& v) const {
        quaternion qConjugate{
            orientation.w,
            -orientation.x,
            -orientation.y,
            -orientation.z
        };

        return rotateVector(v, qConjugate);
    };  
    //------------------------------------------------------------------------------------

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


    //Compute dynamic surface area of the rocket as a function of its orientation
    //This way, drag is realisticlly calculated as the attitude of the rocket changes
    double calculateProjectedArea(const xyz& relativeBody) const {

        double speed = std::sqrt(
            relativeBody.x * relativeBody.x +
            relativeBody.y * relativeBody.y +
            relativeBody.z * relativeBody.z
        );

        if (speed < 1e-8) {
            return frontArea;
        }

        double axialVelocity = std::abs(relativeBody.x);

        double lateralVelocity = std::sqrt(
            relativeBody.y * relativeBody.y + 
            relativeBody.z * relativeBody.z
        );

        double sideArea = length * diameter;

        double projectedArea = frontArea * axialVelocity / speed + sideArea * lateralVelocity / speed;

        return projectedArea;
    }

    //This calculates drag as a vector
    // Calculate aerodynamic drag using orientation-dependent projected area
    aero calculateAerodynamics(double rho) const {

        // Rocket velocity in body coordinates
        xyz velBody = worldToBody(currentVel);

        // Velocity at center of pressure
        xyz rotationVel = crossProduct(angularVelocity, cpOffset);

        xyz cpVelocity{
            velBody.x + rotationVel.x,
            velBody.y + rotationVel.y,
            velBody.z + rotationVel.z
        };

        double lateralSpeed = std::sqrt(
            cpVelocity.y * cpVelocity.y +
            cpVelocity.z * cpVelocity.z
        );

        double axialArea = frontArea;
        double lateralArea = length * diameter;

        // Temporary aerodynamic coefficients
        double CdAxial = 0.8;
        double CdLateral = 1.1;

        // Axial aerodynamic resistance
        double axialForce =
            -0.5 * rho * CdAxial * axialArea
            * std::abs(velBody.x) * velBody.x;

        // Lateral aerodynamic resistance
        xyz lateralForce{
            0.0,
            -0.5 * rho * CdLateral * lateralArea
                * lateralSpeed * cpVelocity.y,

            -0.5 * rho * CdLateral * lateralArea
                * lateralSpeed * cpVelocity.z
        };

        // Total aerodynamic force in body coordinates
        xyz forceBody{
            axialForce,
            lateralForce.y,
            lateralForce.z
        };

        // Aerodynamic torque about center of mass
        xyz torqueBody = crossProduct(cpOffset, lateralForce);

        return {
            bodyToWorld(forceBody),
            torqueBody
        };
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

    //Calculate the orientation of the rocket as a quuaternion through semi-implicit euler integration
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

    //Compute the aerodynamic torque forces acting on the rocket by considering the normal 
    //forces as a function of the rockets orientation
    void calcAeroTorque(const xyz& aeroForce) {

        torque.x = cpOffset.y * aeroForce.z
                    - cpOffset.z * aeroForce.y;

        torque.y = cpOffset.z * aeroForce.x
                    - cpOffset.x * aeroForce.z;

        torque.z = cpOffset.x * aeroForce.y
                    - cpOffset.y * aeroForce.x;
    }

    // Calculates the aerodynamic normal force based on the rocket's
    // orientation, translational velocity, and angular velocity.



    void updateState(const earth& e) {

        angularVel.push_back(angularVelocity);
        orientations.push_back(orientation);

        while (true) {

            //  1) Get gravity and air density as function of altitude
            double gravity = e.getGravity(currentPos.z);
            double rho = e.getAirDensity(currentPos.z);

            //  2) Get aerodynamic forces as function of density (USE 0.8 as Cd for now until we get accurate dynamic estimates)
            aero aerodynamic = calculateAerodynamics(rho);
            torque = aerodynamic.torque;

            //  3) Translational accel
            xyz accel {
                aerodynamic.force.x / mass,
                aerodynamic.force.y / mass,
                aerodynamic.force.z / mass - gravity
            };

            // 4) Rotational Accel
            calculateAngularAccel();


            //  5) Update vel
            currentVel.x += accel.x * dt;
            currentVel.y += accel.y * dt;
            currentVel.z += accel.z * dt;
            updateAngularrVel(dt);

            //  6) Update pos + orientation
            currentPos.x += currentVel.x * dt;
            currentPos.y += currentVel.y * dt;
            currentPos.z += currentVel.z * dt;
            updateOrientation(dt);

            //  7) store updated state
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
    for (size_t i = 0; i<r.pos.size(); i+=100) {
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