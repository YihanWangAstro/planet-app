#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

constexpr bool coll_detect = true;

constexpr double resonance_repeat = 20;

constexpr double sigma_field = 30_kms;

constexpr double sigma_globular = 5_kms;

constexpr double signam_open = 0.1_kms;

constexpr double a_j_field = 1_AU;

constexpr double a_j_globular = 0.1_AU;

constexpr double a_j_open = 1_AU;

constexpr double interact_factor = 0.02;

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

void single_single(ConFile out_file, size_t sim_num, double m_star, double a_j, double sigma, double a_s, double m_s1) {
  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_star);

  double r_s1 = stellar::stellar_radius(stellar::StarType::STAR, m_s1);

  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{m_star, r_d}, jupiter{1_Mj, 1_Rj}, star1{m_s1, r_s1};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    move_to_COM_frame(sun, jupiter);

    auto v_inf = space::random::Maxwellian(sigma);

    auto b_max = scattering::b_max(cluster(sun, jupiter), star1, v_inf, interact_factor);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter), star1, v_inf, b_max, delta);

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
        [&](auto &ptc) { out_file << PACK(i, ptc, jupiter_orbit, in_orbit, b_max, v_inf, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1};

    simulator.run(args);
  }
}

void single_binary(ConFile out_file, size_t sim_num, double m_star, double a_j, double sigma, double a_s, double m_s1) {
  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_star);

  double r_s1 = stellar::stellar_radius(stellar::StarType::STAR, m_s1);

  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{m_star, r_d}, jupiter{1_Mj, 1_Rj}, star1{m_s1, r_s1}, star2{m_s1, r_s1};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    move_to_COM_frame(sun, jupiter);

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, star2);

    auto v_inf = space::random::Maxwellian(sigma);

    auto b_max = scattering::b_max(cluster(sun, jupiter), cluster(star1, star2), v_inf, interact_factor);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter), cluster(star1, star2), v_inf, b_max, delta);

    move_particles(in_orbit, star1, star2);

    move_to_COM_frame(sun, jupiter, star1, star2);

    double end_time = ((E_tot(sun, jupiter, star1, star2) > 0) ? 2.0 : resonance_repeat) *
                      time_to_periapsis(cluster(sun, jupiter), cluster(star1, star2));

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition(collision);
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { out_file << PACK(i, ptc, jupiter_orbit, in_orbit, binary_orbit, b_max, v_inf, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

void binary_single(ConFile out_file, size_t sim_num, double m_star, double a_j, double sigma, double a_s, double m_s1) {
  if (3.7 * a_j > a_s) return;

  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_star);

  double r_s1 = stellar::stellar_radius(stellar::StarType::STAR, m_s1);

  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{m_star, r_d}, sun2{m_star, r_d}, jupiter{1_Mj, 1_Rj}, star1{m_s1, r_s1};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    move_to_COM_frame(sun, jupiter);

    auto binary_orbit = EllipOrbit{sun.mass + jupiter.mass, sun2.mass, a_s, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, sun2);

    move_to_COM_frame(sun, jupiter, sun2);

    auto v_inf = space::random::Maxwellian(sigma);

    auto b_max = scattering::b_max(cluster(sun, jupiter, sun2), star1, v_inf, interact_factor);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter, sun2), star1, v_inf, b_max, delta);

    move_particles(in_orbit, star1);

    move_to_COM_frame(sun, jupiter, star1, sun2);

    double end_time = ((E_tot(sun, jupiter, star1, sun2) > 0) ? 2.0 : resonance_repeat) *
                      time_to_periapsis(cluster(sun, jupiter, sun2), star1);

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition(collision);
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { out_file << PACK(i, ptc, jupiter_orbit, in_orbit, binary_orbit, b_max, v_inf, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, sun2, jupiter, star1};

    simulator.run(args);
  }
}

void binary_binary(ConFile out_file, size_t sim_num, double m_star, double a_j, double sigma, double a_s, double m_s1) {
  if (3.7 * a_j > a_s) return;

  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_star);

  double r_s1 = stellar::stellar_radius(stellar::StarType::STAR, m_s1);

  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{m_star, r_d}, sun2{m_star, r_d}, jupiter{1_Mj, 1_Rj}, star1{m_s1, r_s1}, star2{m_s1, r_s1};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    move_to_COM_frame(sun, jupiter);

    auto binary_orbit = EllipOrbit{sun.mass + jupiter.mass, sun2.mass, a_s, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, sun2);

    move_to_COM_frame(sun, jupiter, sun2);

    auto binary_orbit2 = EllipOrbit{star1.mass, star2.mass, a_s, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(binary_orbit2, star2);

    auto v_inf = space::random::Maxwellian(sigma);

    auto b_max = scattering::b_max(cluster(sun, jupiter, sun2), cluster(star1, star2), v_inf, interact_factor);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter, sun2), cluster(star1, star2), v_inf, b_max, delta);

    move_particles(in_orbit, star1, star2);

    move_to_COM_frame(sun, jupiter, star1, sun2, star2);

    double end_time = ((E_tot(sun, jupiter, star1, sun2, star2) > 0) ? 2.0 : resonance_repeat) *
                      time_to_periapsis(cluster(sun, jupiter, sun2), cluster(star1, star2));

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition(collision);
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc) {
      out_file << PACK(i, ptc, jupiter_orbit, in_orbit, binary_orbit, binary_orbit2, b_max, v_inf, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, sun2, jupiter, star1, star2};

    simulator.run(args);
  }
}

void explore(std::string const &output_name, std::string const &sim_type, size_t sim_num, double m_s, double a_s,
             double m_s1, double aj, double sigma) {
  auto ss_file = make_thread_safe_fstream(output_name + "_" + sim_type + "_ss.txt", std::fstream::out);
  auto sb_file = make_thread_safe_fstream(output_name + "_" + sim_type + "_sb.txt", std::fstream::out);
  auto bs_file = make_thread_safe_fstream(output_name + "_" + sim_type + "_bs.txt", std::fstream::out);
  auto bb_file = make_thread_safe_fstream(output_name + "_" + sim_type + "_bb.txt", std::fstream::out);

  std::vector<std::thread> threads;

  auto thread_num = machine_thread_num / 4;

  for (size_t i = 0; i < thread_num; ++i)
    threads.emplace_back(std::thread(single_single, ss_file, sim_num, m_s, aj, sigma, a_s, m_s1));

  for (size_t i = 0; i < thread_num; ++i)
    threads.emplace_back(std::thread(single_binary, sb_file, sim_num, m_s, aj, sigma, a_s, m_s1));

  for (size_t i = 0; i < thread_num; ++i)
    threads.emplace_back(std::thread(binary_single, bs_file, sim_num, m_s, aj, sigma, a_s, m_s1));

  for (size_t i = 0; i < thread_num; ++i)
    threads.emplace_back(std::thread(binary_binary, bb_file, sim_num, m_s, aj, sigma, a_s, m_s1));

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

  // double const v_factor = sqrt(8 / consts::pi);

  if (sim_type == "open") {
    explore(output_name, "open", sim_num, m_s, a_s, m_s1, a_j_open, signam_open);
  } else if (sim_type == "field") {
    explore(output_name, "filed", sim_num, m_s, a_s, m_s1, a_j_field, sigma_field);
  } else if (sim_type == "globular") {
    explore(output_name, "globular", sim_num, m_s, a_s, m_s1, a_j_globular, sigma_globular);
  } else {
    std::cout << "undefined sim type!\n";
  }
  return 0;
}