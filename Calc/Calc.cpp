#include "Calc.h"

void Calculate::Solve() {
	std::shared_ptr<Parser> p = std::make_shared<Parser>();
	p->read(file);

	Data data(p);

	std::shared_ptr<Solver> solver = FabricSolver::createSolver(data);
	solver->solve();

	Smoothing stress(data, { STRESS });
	stress.solve();

	//Smoothing strain(data, { STRAIN, STRESS });
	//strain.solve();

	VTUWriter writer(data.get_elements(), data.get_nodes());
	writer.write(file + ".vtu");

    try { // change
        std::vector<fs::path> vtu_files;

        // Собираем все VTU файлы в директории
        for (const auto& entry : fs::directory_iterator(std::filesystem::current_path())) {
            if (entry.path().extension() == ".vtu") {
                vtu_files.push_back(entry.path());
            }
        }

        if (vtu_files.empty()) {
            std::cerr << "No VTU files found in directory\n";
            return;
        }

        auto timesteps = parse_timesteps(vtu_files);

        std::sort(timesteps.begin(), timesteps.end(),
                  [](const TimestepFile& a, const TimestepFile& b) {
                      return a.time < b.time;
                  });

        create_pvd_file(timesteps, file + ".pvd");
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return;
    }
}
