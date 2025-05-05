#include "Calc.h"

void Calculate::Solve() {
	std::shared_ptr<Parser> p = std::make_shared<Parser>();
	p->read(file + ".fc");

	Data data(p);

	std::shared_ptr<Solver> solver = FabricSolver::createSolver(data);
	solver->solve();

	Smoothing stress(data, STRESS);
	stress.solve();

	Smoothing strain(data, STRAIN);
	strain.solve();

	VTUWriter writer(data.get_elements(), data.get_nodes());
	writer.write(file + ".vtu");
}
