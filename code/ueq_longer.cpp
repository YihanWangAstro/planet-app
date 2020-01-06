#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

constexpr bool coll_detect = false;

constexpr double resonance_repeat = 10000;

std::array<double, 1> V_INF = {0.1_kms, 1_kms};

std::array<double, 1> AJ = {0.5_AU, 1_AU, 2_AU, 5_AU};

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

void single_single(std::string workdir, size_t sim_num, double m_star, double a_j, double v_inf, double a_s,
                   double m_s1) {
  std::fstream out_file{workdir + ".txt", std::fstream::out};

  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_star);

  double r_s1 = stellar::stellar_radius(stellar::StarType::STAR, m_s1);

  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{m_star, r_d}, jupiter{1_Mj, 1_Rj}, star1{m_s1, r_s1};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    move_to_COM_frame(sun, jupiter);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter), star1, v_inf, delta);

    move_particles(in_orbit, star1);

    move_to_COM_frame(sun, jupiter, star1);

    double end_time =
        ((E_tot(sun, jupiter, star1) > 0) ? 2.0 : resonance_repeat) * time_to_periapsis(cluster(sun, jupiter), star1);

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition(collision);
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { space::display(out_file, i, ptc, jupiter_orbit, in_orbit, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1};

    simulator.run(args);
  }
}

using Task = std::function<void(std::string, size_t, double, double, double, double, double)>;

void explore(Task job, std::string const &output_name, std::string const &sim_type, size_t sim_num, double m_s,
             double a_s, double m_s1) {
  std::vector<std::thread> threads;

  for (auto v : V_INF) {
    for (auto aj : AJ) {
      for (int i = 0; i < 10; ++i) {
        char params[105];
        sprintf(params, "_%.1lf_%.1lf_%d", v / kms, aj, i);
        std::string fname = output_name + "_" + sim_type + params;
        threads.emplace_back(std::thread(job, fname, sim_num, m_s, aj, v, a_s, m_s1));
      }
    }
  }

  for (auto &th : threads) {
    th.join();
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_name;
  std::string sim_type;
  double a_s, m_s, m_s1;

  tools::read_command_line(argc, argv, sim_type, sim_num, output_name, a_s, m_s, m_s1);

  a_s *= unit::AU;
  m_s *= unit::Ms;

  explore(single_single, output_name, "ss", sim_num, m_s, a_s, m_s1);

  return 0;
}