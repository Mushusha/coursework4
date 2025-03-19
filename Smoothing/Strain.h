#pragma once

#include "Smoothing.h"

class Strain : public Smoothing {
public:
    Strain(Data& data, std::vector<double> start,
        std::vector<double> end, int count)
        : Smoothing(data, start, end, count) {}
    Strain(const Strain& other) : Smoothing(other) {} 
    Strain& operator=(const Strain& other) {
        if (this != &other) {
            Smoothing::operator=(other);
        }
        return *this;
    }
    Strain(Strain&& other) noexcept 
        : Smoothing(std::move(other)) {}
    Strain& operator=(Strain&& other) noexcept {
        if (this != &other) {
            Smoothing::operator=(std::move(other));
        }
    }
    ~Strain() = default;

};