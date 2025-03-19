#pragma once

#include "Solver.h"

class Statics : public Solver {
public:
	Statics(Data& data) : Solver(data) {};
    Statics(const Statics& other) : Solver(other) {}
    Statics& operator=(const Statics& other) {
        if (this != &other) {
            Solver::operator=(other);
        }
        return *this;
    }
    Statics(Statics&& other) noexcept
        : Solver(std::move(other)) {}
    Statics& operator=(Statics&& other) noexcept {
        if (this != &other) {
            Solver::operator=(std::move(other));
        }
    }
	virtual ~Statics() = default;

private:
	void calcDisp() final;
};