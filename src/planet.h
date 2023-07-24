#ifndef PWM_MODEL_PLANET_H
#define PWM_MODEL_PLANET_H

#include "json.hpp"
#include <string>
#include <unordered_map>
#include <vector>

namespace PWM{
    namespace Model{
        class planet{
            private:
                //Planet's name
                std::string name;

                //Radius in meters
                double radius;

                //The sidereal otational period in seconds
                double rotPeriod;

                //Average temp in Kelvin
                double aveTemp;

                //Average pressure at sea level in pascal
                double averageSealevelPressure;// 101325 known average for Earth
                
                //Altitude of sea level (0 by convention) in metres
                double referenceHeight;

                //Acceleration due to gravity in the atmospheric environment in metres per second^2, constant for our purposes
                double gravitationalAcceleration;// 9.81 on Earth
                
                //The maximum vapour density (i.e., moisture) in kg per metre^3, used for initialisation.
                double maxVapourDensity;// 0.294271 'max' on Earth, calculated from saturation pressure at 80 deg C
                
                //The maximum wind speed in m per second, used for initialisation.
                double maxWindSpeed;// 40 is an assumption for earth.
                
                //The molar mass of the air, in kg per mole
                double molarMassOfAir;// 0.0289644 kg per mole for earth atmosphere.
                
                //The height of the atmosphere in metres (used in the terrainHeatingEngine)
                double atmosphereHeight;

                //A list of potential terrain types, used for reading in the other maps.
                std::vector<std::string> terrainTypes;

                //The density of the planet terrain, both land and water.
                std::unordered_map<std::string, double> terrainDensity;//On earth, 2600 for soil and 1000 for ocean.

                //The rate at which temperature drops with increasing heigh in kelvin per metre.
                double tempLapseRate;// 0.0098 kelvin per metre of altitude on earth (dry air).

                //The specific heat capacity of the planet's atmosphere.
                double airHeatCapacity; // 1006 joule per kilogram kelvin on earth.

                //The specific heat capacity of the substance that drives the moisture cycle (water on earth)
                double moistureHeatCapacity; //~4000 joule per kilogram kelvin

                //The temperature where fluid turns into gas of the substance that drives the moisture cycle (water on earth)
                double moistureBoilTemp; //373.16 K for water

                //The latent heat of evaporation of moisture in the air
                double moistureLatentHeatofVaporisation; //2257 kJ per kg for water.

                //The specific heat capacity of the planet's terrain, both land and ocean.
                std::unordered_map<std::string, double> terrainSpecHeatCapacity;
                
                //The albedo for the planet's terrain, covering all types.
                std::unordered_map<std::string, double> terrainAlbedo;

                //The emissivity constants for the planet's terrain, differing per type.
                //Used in the atmosphericHeatingEngine
                std::unordered_map<std::string, double> terrainEmissivityConstants;

                //Mapping of if the terrain type is dry (i.e., should take moisture level into account when interacting)
                std::unordered_map<std::string, bool> terrainIsDry;
            public:
                //Default constructor
                planet(double rad = 6371000, double rotP = 86400, double aveTemp = 288.15);

                //Constructor from JSON file
                planet(std::string JSONfile);

                //Write planet data to JSON file.
                int writePlanetFile(std::string file) const;

                //Getters for the various values, will be expanded.
                const std::string getName() const;
                const double getRadius() const;
                const double getRotPeriod() const;
                const double getAveTemp() const;
                const double getAverageSealevelPressure() const;
                const double getReferenceHeight() const;
                const double getGravitationalAcceleration() const;
                const double getMaxVapourDensity() const;
                const double getMaxWindSpeed() const;
                const double getMolarAirMass() const;
                const double getAtmoThickness() const;
                const std::vector<std::string>& getTerrainTypes() const;
                const std::unordered_map<std::string, double>& getTerrainDensity() const;
                const double getTempLapseRate() const;
                const double getAirHeatCapacity() const;
                const double getMoistureHeatCapacity() const;
                const double getMoistureBoilTemp() const;
                const std::unordered_map<std::string, double>& getTerrainHeatCapacity() const;
                const std::unordered_map<std::string, double>& getTerrainAlbedo() const;
                const std::unordered_map<std::string, double>& getTerrainEmissivityConstants() const;
                const std::unordered_map<std::string, bool>& getTerrainDryType() const;
                const double getLatentHeatofVaporisation() const;

                const double getAngularVelocity() const;
                
                //Setters for the various values, not intended to be used.
                void setName(std::string n);
                void setRadius(double rad);
                void setRotPeriod(double rP);
                void setAveTemp(double aTemp);
                void setAverageSealevelPressure(double aslPres);
                void setReferenceHeight(double refHeight);
                void setGravitationalAcceleration(double g);
                void setMaxVapourDensity(double mvd);
                void setMaxWindSpeed(double mws);
                void setMolarAirMass(double mam);
                void setAtmoThickness(double ast);
                void setTerrainTypes(std::vector<std::string> terTypes);
                void setTerrainDensity(std::unordered_map<std::string, double> terDens);
                void setTempLapseRate(double tlr);
                void setAirHeatCap(double ahc);
                void setMoistureHeatCap(double mhc);
                void setMoistureBoilTemp(double mct);
                void setTerrainHeatCapacity(std::unordered_map<std::string, double> terHCP);
                void setTerrainAlbedo(std::unordered_map<std::string, double> terA);
                void setTerrainEmissivityConstants(std::unordered_map<std::string, double> terEC);
                void setTerrainDryType(std::unordered_map<std::string, bool> tdt);
                void setLatentHeatofVaporisation(double lhv);

                bool operator==(const planet& other) const;
                bool operator!=(const planet& other) const;
        };

        //JSON lib interaction functions
        void to_json(nlohmann::json& j, const planet& p);
        void from_json(const nlohmann::json& j, planet& p);
    }
}
#endif //PWM_MODEL_PLANET_H
