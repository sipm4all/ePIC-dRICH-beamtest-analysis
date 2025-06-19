
//  Photon structure
struct Photon
{
    double wavelength; // Wavelength of the photon
    double ref_index;  // Wavelength dependent ref index in lucite
    double x, z;       // Track in the x-z plane only
    double vx, vz;     // Speed components in x and z directions
    bool alive;        // Flag to indicate if the photon is absorbed or not
    bool escaped;      // Flag to indicate if the photon is absorbed or not
};
std::array<double, 2> cherenkov_emission_wavelength = {280, 850};

TH2F *test_hits = new TH2F("test_hits", ";x (mm);z (mm)", 100, -100, 100, 100, -5, 5);
TH2F *test_hits_sensors = new TH2F("test_hits_sensors", ";x (mm);z (mm)", 100, -100, 100, 100, -55, -45);
TH1F *test_hits_sensors_plane = new TH1F("test_hits_sensors_plane", ";x (mm);z (mm)", 260, -175, 175);

//  Simulation parameters
const double PI = 3.141592653589793;
const double n_lucite = 1.49;
const double absorption_length = 1.e4;                     // mm
const double thickness = 2.0;                              // mm
const double radius = 90.0;                                // mm (radial extent)
const double sensor_z = -51.0;                             // mm
const float _c0 = 299792458., speed_of_light = 299792458.; // [m / s]
TF1 *lucite_ref_index = new TF1("lucite_ref_index", "[0] +[3]*TMath::Exp(-[1]*(x-[2]))+[6]*TMath::Exp(-[4]*(x-[5]))", 0, 100);
lucite_ref_index->SetParameters(1.48757, 0.0498538, 275.618, 0.000535379, 0.0103045, -49.8218, 1.65533);

// Utility: generate uniform random double in [0,1)à
inline double rand_uniform()
{
    return static_cast<double>(std::rand()) / RAND_MAX;
}

// Utility: Create a thread-local random engine and normal distribution
inline double rand_gaussian(double mean = 0.0, double stddev = 1.0)
{
    static std::mt19937 generator(std::random_device{}());
    static std::normal_distribution<double> distribution;

    // Update the distribution parameters
    distribution.param(std::normal_distribution<double>::param_type(mean, stddev));
    return distribution(generator);
}

// Generate initial Cherenkov photon in x-z plane
Photon generate_photon(double particle_beta)
{
    // Generates a random wavelength following the spectrum of Cherenkov radiation
    double emission_wavelength = (cherenkov_emission_wavelength[0]) / (1 - (gRandom->Uniform(0, 1)) * (cherenkov_emission_wavelength[1] - cherenkov_emission_wavelength[0]) / (cherenkov_emission_wavelength[1]));

    //  Define the corresponding refractive index associated to the wavelength
    auto current_ref_index = lucite_ref_index->Eval(emission_wavelength);
    bool left_right_orientation = rand_uniform() < 0.5;

    //  Generate the Cherenkov angle
    double cerenkov_angle = std::acos(1.0 / (particle_beta * current_ref_index));

    //  Generate the photon dynamic
    double vx = left_right_orientation ? std::sin(cerenkov_angle) : -std::sin(cerenkov_angle);
    double vz = std::cos(cerenkov_angle);
    return {emission_wavelength, current_ref_index, rand_gaussian(0., 0.1), 2 * rand_uniform(), vx, vz, true};
}

void propagate_photon(Photon &p)
{
    // Propagate the photon in the x-z plane
    auto z_time_step = fabs(p.vz < 0 ? -p.z / p.vz : (thickness - p.z) / p.vz);
    auto x_time_step = fabs(p.vx < 0 ? (radius + p.x) / p.vx : (radius - p.x) / p.vx);
    auto time_step = std::min(1. * z_time_step, 1. * x_time_step);

    // Update position
    auto z_space_step = p.vz * time_step;
    auto x_space_step = p.vx * time_step;
    auto full_space_step = std::sqrt(z_space_step * z_space_step + x_space_step * x_space_step);
    p.z += z_space_step;
    p.x += x_space_step;
    test_hits->Fill(p.x, p.z);

    // Check if the photon is absorbed
    auto absorption_probability = 1.0 - std::exp(-full_space_step / absorption_length);
    if (rand_uniform() < absorption_probability)
    {
        p.alive = false; // Photon is absorbed
        return;
    }

    // Update speed
    //  Check if the photon has a TIR
    double incidence_angle = std::acos(std::abs(p.vz));
    double critical_angle = std::asin(1.0 / p.ref_index);
    if (incidence_angle > critical_angle)
    {
        if (time_step == z_time_step)
            p.vz *= -1; // Reflect off z boundary
        else
        {
            p.vx *= -1; // Reflect off x boundary
            if (rand_uniform() > 0.9)
                p.alive = false; // Photon is absorbed
        }

        // Simulate imperfections with angle perturbation
        auto angle_perturbation = 2 * TMath::Pi() * rand_gaussian(0., 0.05) / 360.0;
        auto new_angle = std::atan2(p.vz, p.vx) + angle_perturbation;
        p.vx = std::cos(new_angle);
        p.vz = std::sin(new_angle);
    }
    else
    {
        p.alive = false; // Photon is absorbed
        if (p.vz < 0)
        {
            p.escaped = true; // Photon escapes
            z_time_step = fabs(p.z + sensor_z / p.vz);
            z_space_step = p.vz * z_time_step;
            x_space_step = p.vx * z_time_step;
            p.x += x_space_step;
            p.z += z_space_step;
            test_hits_sensors->Fill(p.x, p.z);
            test_hits_sensors_plane->Fill(p.x);
        }
    }
}

// Simulate photons through Lucite slice and track those hitting sensor plane
std::vector<std::array<double, 2>> simulate_photons(int num_photons)
{
    std::vector<std::array<double, 2>> hits;
    for (int i = 0; i < num_photons; ++i)
    {
        Photon p = generate_photon(1.);

        while (true)
        {
            propagate_photon(p);

            if (p.escaped)
                hits.push_back({p.x, p.z});

            if (!p.alive)
                break;
        }
    }
    return hits;
}

void test_quick_mc()
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    int num_photons = 100000;
    auto hits = simulate_photons(num_photons);

    std::cout << "Detected " << hits.size() << " photons on the sensor plane (1D projection)." << std::endl;
    for (auto x : hits)
    {
        // std::cout << "x: " << x[0] << " - z: " << x[1] << std::endl;
    }

    TCanvas *c1 = new TCanvas();
    test_hits->Draw("COLZ");
    TCanvas *c2 = new TCanvas();
    test_hits_sensors->Draw("COLZ");

    TCanvas *c3 = new TCanvas();
    test_hits_sensors_plane->Draw();
}