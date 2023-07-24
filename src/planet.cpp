#include <fstream>
#include <iomanip>
#include <iostream>

#include "json.hpp"
#include "planet.h"

namespace PWM{
    namespace Model{
        planet::planet(double rad, double rotP, double aveTemp) : radius(rad), rotPeriod(rotP), aveTemp(aveTemp){

        }

        planet::planet(std::string JSONfile){
            std::ifstream fil(JSONfile);
            if (!fil.good()){
                std::cout << "\033[1;31mError! Planet json file " << JSONfile << " not found!\033[0m" << std::endl;
                return;
            }
            nlohmann::json jsonData;
            try{
                fil >> jsonData;
                *this = jsonData.get<planet>();
            }
            catch (nlohmann::json::exception& e){
                std::cout << "\033[1;31mError! Planet json file not in correct format!\033[0m" << std::endl;
                std::cout << "\033[1;37mPlanet json file not read in.\033[0m" << std::endl;
            }
        }

        int planet::writePlanetFile(std::string file) const{
            std::ofstream fil(file);
            if (!fil.good()){
                std::cerr << "\033[1;31mError! Planet file " << file << " unsuitable for writing!\033[0m" << std::endl;
                return -1;
            }
            try {
                nlohmann::json jsonData = *this;
                fil << std::setw(4) << jsonData << std::endl;
                return 0;
            } catch (nlohmann::json::exception& e) {
                std::cerr << "\033[1;31mError when writing planet json file!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mPlanet json file not written.\033[0m" << std::endl;
                return -2;
            }
        }

        const std::string planet::getName() const{
            return name;
        }

        const double planet::getRadius() const{
            return radius;
        }

        const double planet::getRotPeriod() const{
            return rotPeriod;
        }

        const double planet::getAveTemp() const{
            return aveTemp;
        }

        const double planet::getAverageSealevelPressure() const{
            return averageSealevelPressure;
        }

        const double planet::getReferenceHeight() const{
            return referenceHeight;
        }

        const double planet::getGravitationalAcceleration() const{
            return gravitationalAcceleration;
        }

        const double planet::getMaxVapourDensity() const{
            return maxVapourDensity;
        }

        const double planet::getMaxWindSpeed() const{
            return maxWindSpeed;
        }

        const double planet::getMolarAirMass() const{
            return molarMassOfAir;
        }

        const double planet::getAtmoThickness() const{
            return atmosphereHeight;
        }

        const std::vector<std::string>& planet::getTerrainTypes() const{
            return terrainTypes;
        }

        const std::unordered_map<std::string, double>& planet::getTerrainDensity() const{
            return terrainDensity;
        }

        const double planet::getTempLapseRate() const{
            return tempLapseRate;
        }

        const double planet::getAirHeatCapacity() const{
            return airHeatCapacity;
        }

        const double planet::getMoistureHeatCapacity() const{
            return moistureHeatCapacity;
        }

        const double planet::getMoistureBoilTemp() const{
            return moistureBoilTemp;
        }

        const std::unordered_map<std::string, double>& planet::getTerrainHeatCapacity() const{
            return terrainSpecHeatCapacity;
        }

        const std::unordered_map<std::string, double>& planet::getTerrainAlbedo() const{
            return terrainAlbedo;
        }

        const std::unordered_map<std::string, double>& planet::getTerrainEmissivityConstants() const{
            return terrainEmissivityConstants;
        }

        const std::unordered_map<std::string, bool>& planet::getTerrainDryType() const{
            return terrainIsDry;
        }

        const double planet::getLatentHeatofVaporisation() const{
            return moistureLatentHeatofVaporisation;
        }

        const double planet::getAngularVelocity() const{
            #if __cplusplus > 201703L
                return 2 * std::numbers::pi / rotPeriod;
            #else
                return 2 * 3.1415926535 / rotPeriod;
            #endif
        }

        void planet::setName(std::string n){
            name = n;
        }

        void planet::setRadius(double rad){
            radius = rad;
        }

        void planet::setRotPeriod(double rP){
            rotPeriod = rP;
        }
        
        void planet::setAveTemp(double aTemp){
            aveTemp = aTemp;
        }

        void planet::setAverageSealevelPressure(double aslPres){
            averageSealevelPressure = aslPres;
        }

        void planet::setReferenceHeight(double refHeight){
            referenceHeight = refHeight;
        }

        void planet::setGravitationalAcceleration(double g){
            gravitationalAcceleration = g;
        }

        void planet::setMaxVapourDensity(double mvd){
            maxVapourDensity = mvd;
        }

        void planet::setMaxWindSpeed(double mws){
            maxWindSpeed = mws;
        }

        void planet::setMolarAirMass(double mam){
            molarMassOfAir = mam;
        }

        void planet::setAtmoThickness(double ast){
            atmosphereHeight = ast;
        }

        void planet::setTerrainTypes(std::vector<std::string> terTypes){
            terrainTypes = terTypes;
        }

        void planet::setTerrainDensity(std::unordered_map<std::string, double> terDens){
            terrainDensity = terDens;
        }

        void planet::setTempLapseRate(double tlr){
            tempLapseRate = tlr;
        }

        void planet::setAirHeatCap(double ahc){
            airHeatCapacity = ahc;
        }

        void planet::setMoistureHeatCap(double mhc){
            moistureHeatCapacity = mhc;
        }

        void planet::setMoistureBoilTemp(double mct)        {
            moistureBoilTemp = mct;
        }

        void planet::setTerrainHeatCapacity(std::unordered_map<std::string, double> terHCP){
            terrainSpecHeatCapacity = terHCP;
        }

        void planet::setTerrainAlbedo(std::unordered_map<std::string, double> terA){
            terrainAlbedo = terA;
        }

        void planet::setTerrainEmissivityConstants(std::unordered_map<std::string, double> terEC){
            terrainEmissivityConstants = terEC;
        }

        void planet::setTerrainDryType(std::unordered_map<std::string, bool> tdt){
            terrainIsDry = tdt;
        }

        void planet::setLatentHeatofVaporisation(double lhv){
            moistureLatentHeatofVaporisation = lhv;
        }

        bool planet::operator==(const planet& other) const{
            if (this->getRadius() != other.getRadius())
                return false;
            if (this->getRotPeriod() != other.getRotPeriod())
                return false;
            if (this->getAveTemp() != other.getAveTemp())
                return false;
            if (this->getAverageSealevelPressure() != other.getAverageSealevelPressure())
                return false;
            if (this->getReferenceHeight() != other.getReferenceHeight())
                return false;
            if (this->getGravitationalAcceleration() != other.getGravitationalAcceleration())
                return false;
            if (this->getMaxVapourDensity() != other.getMaxVapourDensity())
                return false;
            if (this->getMaxWindSpeed() != other.getMaxWindSpeed())
                return false;
            if (this->getMolarAirMass() != other.getMolarAirMass())
                return false;
            if (this->getAtmoThickness() != other.getAtmoThickness())
                return false;
            if (this->getTerrainTypes() != other.getTerrainTypes())
                return false;
            if (this->getTerrainDensity() != other.getTerrainDensity())
                return false;
            if (this->getTempLapseRate() != other.getTempLapseRate())
                return false;
            if (this->getAirHeatCapacity() != other.getAirHeatCapacity())
                return false;
            if (this->getMoistureHeatCapacity() != other.getMoistureHeatCapacity())
                return false;
            if (this->getMoistureBoilTemp() != other.getMoistureBoilTemp())
                return false;
            if (this->getTerrainHeatCapacity() != other.getTerrainHeatCapacity())
                return false;
            if (this->getTerrainAlbedo() != other.getTerrainAlbedo())
                return false;
            if (this->getTerrainEmissivityConstants() != other.getTerrainEmissivityConstants())
                return false;
            if (this->getTerrainDryType() != other.getTerrainDryType())
                return false;
            if (this->getLatentHeatofVaporisation() != other.getLatentHeatofVaporisation())
                return false;
            return true;
        }

        bool planet::operator!=(const planet& other) const{
            return !(*this == other);
        }

        void to_json(nlohmann::json& j, const planet& p){
            j = nlohmann::json{{"name", p.getName()},
                               {"radius", p.getRadius()},
                               {"rotPeriod", p.getRotPeriod()},
                               {"aveTemp", p.getAveTemp()},
                               {"averageSeaLevelPressure", p.getAverageSealevelPressure()},
                               {"referenceHeight", p.getReferenceHeight()},
                               {"gravitationalAcceleration", p.getGravitationalAcceleration()},
                               {"maxVapourDensity", p.getMaxVapourDensity()},
                               {"maxWindSpeed", p.getMaxWindSpeed()},
                               {"molarAirMass", p.getMolarAirMass()},
                               {"atmoThickness", p.getAtmoThickness()},
                               {"terrainTypes", p.getTerrainTypes()},
                               {"terrainDensity", p.getTerrainDensity()},
                               {"tempLapseRate", p.getTempLapseRate()},
                               {"airSpecHeatCapacity", p.getAirHeatCapacity()},
                               {"moistureSpecHeatCapacity", p.getMoistureHeatCapacity()},
                               {"moistureBoilTemp", p.getMoistureBoilTemp()},
                               {"terrainSpecHeatCapacity", p.getTerrainHeatCapacity()},
                               {"terrainAlbedo", p.getTerrainAlbedo()},
                               {"terrainEmissivityConstants", p.getTerrainEmissivityConstants()},
                               {"terrainIsDry", p.getTerrainDryType()},
                               {"latentHeatofVaporisation", p.getLatentHeatofVaporisation()}};
        }

        void from_json(const nlohmann::json& j, planet& p){
            p.setName(j.at("name").get<std::string>());
            p.setRadius(j.at("radius").get<double>());
            p.setRotPeriod(j.at("rotPeriod").get<double>());
            p.setAveTemp(j.at("aveTemp").get<double>());
            p.setAverageSealevelPressure(j.at("averageSeaLevelPressure").get<double>());
            p.setReferenceHeight(j.at("referenceHeight").get<double>());
            p.setGravitationalAcceleration(j.at("gravitationalAcceleration").get<double>());
            p.setMaxVapourDensity(j.at("maxVapourDensity").get<double>());
            p.setMaxWindSpeed(j.at("maxWindSpeed").get<double>());
            p.setMolarAirMass(j.at("molarAirMass").get<double>());
            p.setAtmoThickness(j.at("atmoThickness").get<double>());
            p.setTerrainTypes(j.at("terrainTypes").get<std::vector<std::string>>());
            p.setTerrainDensity(j.at("terrainDensity").get<std::unordered_map<std::string, double>>());
            p.setTempLapseRate(j.at("tempLapseRate").get<double>());
            p.setAirHeatCap(j.at("airSpecHeatCapacity").get<double>());
            p.setMoistureHeatCap(j.at("moistureSpecHeatCapacity").get<double>());
            p.setMoistureBoilTemp(j.at("moistureBoilTemp").get<double>());
            p.setTerrainHeatCapacity(j.at("terrainSpecHeatCapacity").get<std::unordered_map<std::string, double>>());
            p.setTerrainAlbedo(j.at("terrainAlbedo").get<std::unordered_map<std::string, double>>());
            p.setTerrainEmissivityConstants(j.at("terrainEmissivityConstants").get<std::unordered_map<std::string, double>>());
            p.setTerrainDryType(j.at("terrainIsDry").get<std::unordered_map<std::string, bool>>());
            p.setLatentHeatofVaporisation(j.at("latentHeatofVaporisation").get<double>());
        }
    }
}
