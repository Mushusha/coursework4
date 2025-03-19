#pragma once

#include "Smoothing.h"

class Stress : public Smoothing {
public:
    Stress(Data& data, std::vector<double> start,
        std::vector<double> end, int count)
        : Smoothing(data, start, end, count) {}
    Stress(const Stress& other) : Smoothing(other) {}
    Stress& operator=(const Stress& other) {
        if (this != &other) {
            Smoothing::operator=(other);
        }
        return *this;
    }
    Stress(Stress&& other) noexcept
        : Smoothing(std::move(other)) {}
    Stress& operator=(Stress&& other) noexcept {
        if (this != &other) {
            Smoothing::operator=(std::move(other));
        }
    }
    ~Stress() = default;

};