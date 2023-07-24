#ifndef PWM_UTILS_ATMO
#define PWM_UTILS_ATMO

#define universalGasConstant 8.3144598 //Joule per (mole·kelvin)
#define stefanBoltzmannConstant 0.0000000567037 // Watt per (meter^2 Kelvin^4)

#include <cmath>
#include "planet.h"

namespace PWM{
    namespace Utils{
        /**
         * Calculates pressure at a given height with the barometric formula.
         * See https://en.wikipedia.org/wiki/Barometric_formula.
         */
        template<typename T>
    inline T altitudeAdjustedPressure(T& height, const std::shared_ptr<PWM::Model::planet> P){
            return P->getAverageSealevelPressure() * std::pow(((P->getAveTemp() + P->getTempLapseRate() * (height - P->getReferenceHeight())) / P->getAveTemp()), ((-P->getGravitationalAcceleration() * P->getMolarAirMass()) / (universalGasConstant * P->getTempLapseRate())));
        }

        /**
         * Calculate the density of air at a given height with barometric formula.
         * See https://en.wikipedia.org/wiki/Barometric_formula#Density_equations.
         */
        template<typename T>
    inline  T altitudeAdjustedDensity(T& height, const std::shared_ptr<PWM::Model::planet> P){
            return (altitudeAdjustedPressure(height, P) * P->getMolarAirMass()) / (universalGasConstant * P->getAveTemp() * (1 - ((P->getTempLapseRate() * (height - P->getReferenceHeight())) / P->getAveTemp())));
        }

        /**
         * Calculate the saturation vapour pressure for a given temperature.
         * Using the Arden Buck equations (https:://en.wikipedia.org/wiki/Arden_Buck_equation)
         * since Tetens is not sufficiently accurate below 0.
         * Arden-Buck optimised for temps between 193 and 323 K.
         */
        template<typename T>
    inline T saturationVapourPressure(T& tempinKelvin){
            T tempinCelsius = tempinKelvin - 273.15;
            T ps = 0;
            if (tempinCelsius < 0)
                ps = 611.15 * std::exp((23.036 - tempinCelsius / 333.7) * (tempinCelsius / (279.82 + tempinCelsius)));
            else
                ps = 611.21 * std::exp((18.678 - tempinCelsius / 234.5) * (tempinCelsius / (257.14 + tempinCelsius)));
            return ps;
        }

        /**
         * Calculate the density of air at a given pressure and temperature (ideal gas law).
         */
        template<typename T>
    inline T pressureAdjustedDensity(T& pres, T& temp, const std::shared_ptr<PWM::Model::planet> P){
            return pres / ((universalGasConstant / P->getMolarAirMass()) * temp);
        }
    }
}
#endif
