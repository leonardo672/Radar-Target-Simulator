#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;

double c = 300000000; // Speed of light in m/s
double lambda = 0.055; // Wavelength in meters
double f0 = c / lambda; // Frequency in Hz

// double Rj = sqrt(element.delta_x_j * element.delta_x_j + element.delta_y_j * element.delta_y_j + element.delta_z_j * element.delta_z_j);

double delta_fe = 2 * M_PI; // Assume this value is 2*pi

struct TargetElement {
    double x_i, y_i, z_i;
    double x_j, y_j, z_j;
    double delta_x_j, delta_y_j, delta_z_j;
};

struct TargetParameters {
    std::string purpose;
    double max_size;
    int num_elements;
    double Maximum_target_length;
    double speed;
    std::vector<TargetElement> elements;
    double target_range; // Added target range attribute1
    double Delta_Target_Range;
};

// Function to calculate target range
double calculateElementRange(const TargetElement& element) {
    return sqrt(element.x_i * element.x_i + element.y_i * element.y_i + element.z_i * element.z_i);
}

double calculateTargetRange(const TargetParameters& target) {
    double maxRange = 0.0;

    for (const auto& element : target.elements) {
        double elementRange = calculateElementRange(element);
        if (elementRange > maxRange) {
            maxRange = elementRange;
        }
    }

    return maxRange;
}

double calculateDeltaElementRange(const TargetElement& element) {
    return sqrt(element.delta_x_j * element.delta_x_j + element.delta_y_j * element.delta_y_j + element.delta_z_j * element.delta_z_j);
}

double calculateDeltaTargetRange(const TargetParameters& target) {
    double maxDeltaRange = 0.0;

    for (const auto& element : target.elements) {
        double deltaElementRange = calculateDeltaElementRange(element);
        if (deltaElementRange > maxDeltaRange) {
            maxDeltaRange = deltaElementRange;
        }
    }

    return maxDeltaRange;
}

void inputTargetParameters(TargetParameters& target) {

    //std::cin >> target.purpose;

    //std::cout << "Maximum size of the target (m): ";
    //std::cin >> target.max_size;
    std::cout << "Maximum target length for all observation angles: " << endl;
    cin >> target.Maximum_target_length;

    std::cout << "Required number of target resolution elements to ensure a given probability classification: ";
    std::cin >> target.num_elements;

    std::cout << "Target speed (m/s): ";
    std::cin >> target.speed;

    target.elements.resize(target.num_elements);

    int n = 1;
    for (int i = 0; i < n; ++i) {
        std::cout << "Coordinate of the point1-th element on the X-axis (m): ";
        std::cin >> target.elements[i].x_i;

        std::cout << "Coordinate of the point1-th element on the Y-axis (m): ";
        std::cin >> target.elements[i].y_i;

        std::cout << "Coordinate of the point1-th element on the Z-axis (m): ";
        std::cin >> target.elements[i].z_i;

        std::cout << "Coordinate of the point2-th element (j-th) on the X-axis (m): ";
        std::cin >> target.elements[i].x_j;

        std::cout << "Coordinate of the point2-th element (j-th) on the Y-axis (m): ";
        std::cin >> target.elements[i].y_j;

        std::cout << "Coordinate of the point2-th element (j-th) on the Z-axis (m): ";
        std::cin >> target.elements[i].z_j;

        std::cout << "Coordinate difference of the point3-th element (j-th) on the X-axis (m): ";
        std::cin >> target.elements[i].delta_x_j;

        std::cout << "Coordinate difference of the point3-th element (j-th) on the Y-axis (m): ";
        std::cin >> target.elements[i].delta_y_j;

        std::cout << "Coordinate difference of the point3-th element (j-th) on the Z-axis (m): ";
        std::cin >> target.elements[i].delta_z_j;
    }
}

struct RadarParameters {
    double initial_phase;
    double angle_velocity_sighting;
    double speed_of_light;
    double amplitude_emitted_signal;
    double amplitude_echo_signal;
    double modulation_period;
    double wavelength;
    double phase_difference_elements;
    double initial_frequency;
};

void inputRadarParameters(RadarParameters& radar) {
    cout << "Radar Station Parameters:" << endl;
    std::cout << "Initial phase (degrees): ";
    std::cin >> radar.initial_phase;

    std::cout << "Angle between velocity and sighting (degrees): ";
    std::cin >> radar.angle_velocity_sighting;

    std::cout << "Speed of light (m/s): ";
    std::cin >> radar.speed_of_light;

    std::cout << "Amplitude of the emitted signal: ";
    std::cin >> radar.amplitude_emitted_signal;

    std::cout << "Amplitude of the echo signal: ";
    std::cin >> radar.amplitude_echo_signal;

    std::cout << "Modulation period (seconds): ";
    std::cin >> radar.modulation_period;

    std::cout << "Wavelength (m): ";
    std::cin >> radar.wavelength;

    std::cout << "Phase difference between elements: ";
    std::cin >> radar.phase_difference_elements;

    // Calculate initial frequency based on provided rules (example)
  //  const double base_frequency = 5.0; // Base frequency in GHz (example)
   // const double conversion_factor = 1e9; // GHz to Hz conversion factor
    //double c = 300000000;
    radar.initial_frequency = radar.speed_of_light / radar.wavelength;

    std::cout << "Calculated initial frequency: " << radar.initial_frequency << " Hz" << std::endl; // Output calculated frequency
}

double calculateDeltaDx(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {
   
    //double delta_x_j
   // double f0 = radar.speed_of_light / lambda;
   // double delta_dx = pow(element.x_i + element.delta_x_j, 2) - 4 * ((c * calculateElementRange(element) * delta_fe) / (M_PI * f0));
    double delta_dx = pow(target.elements[0].x_i + target.elements[0].delta_x_j, 2) - 4 * ((c * calculateTargetRange(target) * delta_fe) / (M_PI * f0));
    return delta_dx;
}

double d1_x(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double d1_x_c = ((target.elements[0].x_i + target.elements[0].delta_x_j) - sqrt(calculateDeltaDx(target, element, radar))) / 2;
    return d1_x_c;
}

double d2_x(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double d2_x_c = ((target.elements[0].x_i + target.elements[0].delta_x_j) + sqrt(calculateDeltaDx(target, element, radar))) / 2;
    return d2_x_c;
}

double calculateDeltaDZ(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double delta_dz = pow((target.elements[0].z_i + target.elements[0].delta_z_j), 2) - 4 * (c * calculateTargetRange(target) * delta_fe) / (M_PI * f0);
    return delta_dz;
}

double d1_Z(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double d1_z_c = ((target.elements[0].z_i + target.elements[0].delta_z_j) - sqrt(calculateDeltaDZ(target, element, radar))) / 2;
    return d1_z_c;
}

double d2_Z(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double d2_z_c = ((target.elements[0].z_i + target.elements[0].delta_z_j) + sqrt(calculateDeltaDZ(target, element, radar))) / 2;
    return d2_z_c;
}

double delta_fc(const TargetParameters& target) {

    double delta_fc_c = (c * target.num_elements) / (2 * target.Maximum_target_length);
    return delta_fc_c;
}

double d_r(const TargetParameters& target) {

    double d_r_c = c / (2 * delta_fc(target));
    return d_r_c;
}

double T_c(const TargetParameters& target, RadarParameters& radar) {

    double T_c_c = (calculateTargetRange(target) * lambda * target.num_elements) / (2 * sin(radar.angle_velocity_sighting) * target.Maximum_target_length);
    return T_c_c;
}

double R1_j(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double R1_j_c = sqrt(((target.elements[0].x_i + target.elements[0].delta_x_j), 2) + pow(target.elements[0].y_i, 2) + pow(target.elements[0].z_i, 2) + pow(target.elements[0].x_i - target.elements[0].x_j, 2) + pow(target.elements[0].y_i - target.elements[0].y_j, 2) + pow(target.elements[0].z_i - target.elements[0].z_j, 2));
        
    return R1_j_c;
}

double R2_j(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {

    double R2_j_c = sqrt(pow(target.elements[0].x_i + target.elements[0].delta_x_j - pow(d1_x(target, element, radar), 2) , 2) + pow(target.elements[0].y_i, 2) + pow(target.elements[0].z_i, 2) + pow(target.elements[0].x_i - target.elements[0].x_j, 2) + pow(target.elements[0].y_i - target.elements[0].y_j, 2) + pow(target.elements[0].z_i - target.elements[0].z_j, 2));

    return R2_j_c;
}

double deltafe_112j(const TargetParameters& target, const TargetElement& element, RadarParameters& radar) {
    double d1_x_c = d1_x(target, target.elements[0], radar);
    double R1_j_C = R1_j(target, target.elements[0], radar);
    double R2_j_C = R2_j(target, target.elements[0], radar);

    //double R2_j_c = ((8 * M_PI * f0 * d1_x(target, element, radar) * (target.elements[0].x_i + target.elements[0].delta_x_j)) - (pow(d1_x(target, element, radar), 2)) / (c * R1_j(target, element, radar) + R2_j(target, element, radar)));
   // double R2_j_c = (8 * M_PI * f0 * d1_x(target, element, radar) * (target.elements[0].x_i + target.elements[0].delta_x_j) - pow(d1_x(target, element, radar), 2)) / (c * R1_j(target, element, radar) + R2_j(target, element, radar));
  //    double R2_j_c = 8 * M_PI * f0 * d1_x(target, element, radar) * target.elements[0].x_i + target.elements[0].delta_x_j - pow(d1_x(target, element, radar), 2) / c * R1_j(target, element, radar) + R2_j(target, element, radar);
       // double R2_j_c = (8 * M_PI * f0 * d1_x(target, element, radar) * (target.elements[0].x_i + target.elements[0].delta_x_j) - pow(d1_x(target, element, radar), 2) / (c * R1_j(target, element, radar) + R2_j(target, element, radar)));
        double deltafe_112j = ((8 * M_PI * f0 * d1_x_c * (target.elements[0].x_i + target.elements[0].delta_x_j)) - pow(d1_x_c, 2)) / (c * (R1_j_C + R2_j_C));
        return deltafe_112j;
}

/*double calculateDeltaDx(const TargetElement& element, double initial_frequency) {
    double c = 300000000; // Speed of light in m/s
    double delta_fe = 2 * M_PI;
    double lambda = 0.055;
    double f0 = c / lambda;
    double delta_dx = pow(element.x_i + element.delta_x_j, 2) - 4 * ((c * calculateElementRange(element) * delta_fe) / (M_PI * f0));
    return delta_dx;
}*/

int Number_of_Targets(int NoTargets) {

    for (int i = 1; i <= NoTargets; i++)
    {
        std::cout << "Target " << i << ":" << endl;
        TargetParameters target;
        inputTargetParameters(target);
        // Calculate the target range based on its elements
        target.target_range = calculateTargetRange(target);
        // Use the collected target parameters as needed
        std::cout << "Target Range: " << target.target_range << " meters" << std::endl;

        target.Delta_Target_Range = calculateDeltaTargetRange(target);
        std::cout << "Delta Target Range: " << target.Delta_Target_Range << " meters" << std::endl;

        TargetElement element;
        RadarParameters radar;
        inputRadarParameters(radar);
        double initial_frequency = radar.speed_of_light / radar.wavelength;
        if (!target.elements.empty()) {
            double deltaDx = calculateDeltaDx(target, target.elements[0], radar);
            cout << "Delta_dx, X-axis: " << deltaDx << endl;
        }
        double d1_x_c = d1_x(target, target.elements[0], radar);
        cout << "Maximum distance d between X-axis antenna element = " << d1_x_c << endl;
        //  target.target_range = calculateElementRange(target.elements[0]); // Using the first element's coordinates
        //  std::cout << "Target Range: " << target.target_range << " meters" << std::endl;
        double d2_x_c = d2_x(target, target.elements[0], radar);
        cout << "Minimum distance d between X-axis antenna element = " << d2_x_c << endl;

        double deltaDz = calculateDeltaDZ(target, target.elements[0], radar);
        cout << "Delta_dZ, Z-axis: " << deltaDz << endl;

        double d1_z_c = d1_Z(target, target.elements[0], radar);
        cout << "Maximum distance d between Z-axis antenna element = " << d1_z_c << endl;
        //  target.target_range = calculateElementRange(target.elements[0]); // Using the first element's coordinates
        //  std::cout << "Target Range: " << target.target_range << " meters" << std::endl;
        double d2_z_c = d2_Z(target, target.elements[0], radar);
        cout << "Minimum distance d between Z-axis antenna element = " << d2_z_c << endl;

        double delta_fc_C = delta_fc(target);
        cout << "Width of the spectrum of the probing signal = " << delta_fc_C << endl;

        double d_r_C = d_r(target);
        cout << "Radar resolution by range = " << d_r_C << endl;

        double T_c_C = T_c(target, radar);
        cout << "Antenna aperture synthesis time sec = " << T_c_C << endl;

        cout << "expressions for the phase difference for the radar image elements of the j-th point reflector " << endl
            << "of the first and second, as well as the first and third elements of the antenna array:" << endl;
        double R1_j_C = R1_j(target, target.elements[0], radar);
        cout << fixed << setprecision(5) << "R1_j = " << R1_j_C << endl;

        double R2_j_C = R2_j(target, target.elements[0], radar);
        cout << fixed << setprecision(5) << "R2_j = " << R2_j_C << endl;

        double deltafe_112j_C = deltafe_112j(target, target.elements[0], radar);
        cout << fixed << setprecision(5) << "deltafe_112j = " << deltafe_112j_C << endl;
    }

    return NoTargets;
}

int main() {


    int NoTargets;
    cout << "The number of the targets:" << endl;
    cin >> NoTargets;
    Number_of_Targets(NoTargets);



    return 0;
}