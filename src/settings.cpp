#include <iomanip>
#include "settings.h"

namespace PWM{
    namespace Utils{
        settings::settings(){
        }

        settings::settings(std::string file){
            std::ifstream fil(file);
            if (!fil.good()){
                std::cerr << "\033[1;31mError! Settings file " << file << " not found!\033[0m" << std::endl;
                return;
            }
            nlohmann::json jsonData;
            try{
                fil >> jsonData;
                *this = jsonData.get<settings>();
            }
            catch (nlohmann::json::exception& e){
                std::cerr << "\033[1;31mError! Settings json file not in correct format!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mSettings json file not read in.\033[0m" << std::endl;
            }
        }

        int settings::writeSettings(std::string file) const{
            std::ofstream fil(file);
            if (!fil.good()){
                std::cerr << "\033[1;31mError! Settings file " << file << " unsuitable for writing!\033[0m" << std::endl;
                return -1;
            }
            try {
                nlohmann::json jsonData = *this;
                fil << std::setw(4) << jsonData << std::endl;
                return 0;
            } catch (nlohmann::json::exception& e) {
                std::cerr << "\033[1;31mError when writing settings json file!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mSettings json file not written.\033[0m" << std::endl;
                return -2;
            }
        }

        bool settings::operator==(const settings& other) const{
            if (maxLayers != other.maxLayers)
                return false;
            if (coriolis != other.coriolis)
                return false;
            if (aheSurfaceFluxRatio != other.aheSurfaceFluxRatio)
                return false;
            if (aheThreshold != other.aheThreshold)
                return false;
            if (hdeCoefficient != other.hdeCoefficient)
                return false;
            if (pteDoEvaporation != other.pteDoEvaporation)
                return false;
            if (pteDoCondensation != other.pteDoCondensation)
                return false;
            if (pteDoEvaporationHeatExchange != other.pteDoEvaporationHeatExchange)
                return false;
            if (pteDoCondensationHeatExchange != other.pteDoCondensationHeatExchange)
                return false;
            if (pteDoParticleCondensation != other.pteDoParticleCondensation)
                return false;
            if (pteEvaporationCoefficient != other.pteEvaporationCoefficient)
                return false;
            if (pteCondensationCoefficient != other.pteCondensationCoefficient)
                return false;
            if (theHeatLossCoefficient != other.theHeatLossCoefficient)
                return false;
            if (vceConvectCouplingCoeff != other.vceConvectCouplingCoeff)
                return false;
            if (vcePressureCouplingCoeff != other.vcePressureCouplingCoeff)
                return false;
            if (vcePressureUpliftCoeff != other.vcePressureUpliftCoeff)
                return false;
            if (vcePerformPressureUplift != other.vcePerformPressureUplift)
            if (weDt != other.weDt)
                return false;
            if (P != other.P)
                return false;
            if (suns != other.suns)
                return false;
            if (wmThreshold != other.wmThreshold)
                return false;
            if (atmoTopHeight != other.atmoTopHeight)
                return false;
            if (actualSize != other.actualSize)
                return false;
            if (terrainCellWidth != other.terrainCellWidth)
                return false;
            if (plumeDT != other.plumeDT)
                return false;
            if (subspheres != other.subspheres)
                return false;
            if (subsubspheres != other.subsubspheres)
                return false;
            if (initialPlumeSpeed != other.initialPlumeSpeed)
                return false;
            if (initialPlumeDensity != other.initialPlumeDensity)
                return false;
            if (ventRay != other.ventRay)
                return false;
            if (ventAltitude != other.ventAltitude)
                return false;
            if (initialEjectionAngle != other.initialEjectionAngle)
                return false;

            return true;
        }

        bool settings::operator!=(const settings& other) const{
            return !(*this == other);
        }

        void to_json(nlohmann::json& j, const settings& s){
            j = nlohmann::json{
                {"coriolis", s.coriolis},
                {"aheSurfaceFluxRatio", s.aheSurfaceFluxRatio},
                {"aheThreshold", s.aheThreshold},
                {"hdeCoefficient", s.hdeCoefficient},
                {"pteDoEvaporation", s.pteDoEvaporation},
                {"pteDoCondensation", s.pteDoCondensation},
                {"pteDoEvaporationHeatExchange", s.pteDoEvaporationHeatExchange},
                {"pteDoCondensationHeatExchange", s.pteDoCondensationHeatExchange},
                {"pteDoParticleCondensation", s.pteDoParticleCondensation},
                {"pteEvaporationCoefficient", s.pteEvaporationCoefficient},
                {"pteCondensationCoefficient", s.pteCondensationCoefficient},
                {"theHeatLossCoefficient", s.theHeatLossCoefficient},
                {"vceConvectCouplingCoeff", s.vceConvectCouplingCoeff},
                {"vcePressureCouplingCoeff", s.vcePressureCouplingCoeff},
                {"vcePressureUpliftCoeff", s.vcePressureUpliftCoeff},
                {"vcePerformPressureUplift", s.vcePerformPressureUplift},
                {"weDt", s.weDt},
                {"P", s.P},
                {"suns", s.suns},
                {"wmThreshold", s.wmThreshold},
                {"atmoTopHeight", s.atmoTopHeight},
                {"actualSize", s.actualSize},
                {"maxLayers", s.maxLayers},
                {"gridWidth", s.gridWidth},
                {"terrainCellWidth", s.terrainCellWidth},
                {"terrainFile", s.terrainFile},
                {"smokeTextureFile", s.smokeTextureFile},
                {"plumeDT", s.plumeDT},
                {"subspheres", s.subspheres},
                {"subsubspheres", s.subsubspheres},
                {"initialPlumeSpeed", s.initialPlumeSpeed},
                {"initialPlumeDensity", s.initialPlumeDensity},
                {"ventRay", s.ventRay},
                {"ventAltitude", s.ventAltitude},
                {"initialEjectionAngle", s.initialEjectionAngle}
            };
        }

        void from_json(const nlohmann::json& j, settings& s){
            s.coriolis = j.at("coriolis").get<bool>();
            s.aheSurfaceFluxRatio = j.at("aheSurfaceFluxRatio").get<float>();
            s.aheThreshold = j.at("aheThreshold").get<float>();
            s.hdeCoefficient = j.at("hdeCoefficient").get<float>();
            s.pteDoEvaporation = j.at("pteDoEvaporation").get<bool>();
            s.pteDoCondensation = j.at("pteDoCondensation").get<bool>();
            s.pteDoEvaporationHeatExchange = j.at("pteDoEvaporationHeatExchange").get<bool>();
            s.pteDoCondensationHeatExchange = j.at("pteDoCondensationHeatExchange").get<bool>();
            s.pteDoParticleCondensation = j.at("pteDoParticleCondensation").get<bool>();
            s.pteEvaporationCoefficient = j.at("pteEvaporationCoefficient").get<float>();
            s.pteCondensationCoefficient = j.at("pteCondensationCoefficient").get<float>();
            s.theHeatLossCoefficient = j.at("theHeatLossCoefficient").get<float>();
            s.vceConvectCouplingCoeff = j.at("vceConvectCouplingCoeff").get<float>();
            s.vcePressureCouplingCoeff = j.at("vcePressureCouplingCoeff").get<float>();
            s.vcePressureUpliftCoeff = j.at("vcePressureUpliftCoeff").get<float>();
            s.vcePerformPressureUplift = j.at("vcePerformPressureUplift").get<bool>();
            s.weDt = j.at("weDt").get<float>();
            s.P = j.at("P").get<PWM::Model::planet>();
            s.suns = j.at("suns").get<std::vector<PWM::Model::sun<float>>>();
            s.wmThreshold = j.at("wmThreshold").get<float>();
            s.atmoTopHeight = j.at("atmoTopHeight").get<float>();
            s.actualSize = j.at("actualSize").get<float>();
            s.gridWidth = j.at("gridWidth").get<int>();
            s.terrainCellWidth = j.at("terrainCellWidth").get<int>();
            s.maxLayers = j.at("maxLayers").get<int>();
            s.terrainFile = j.at("terrainFile").get<std::string>();
            s.smokeTextureFile = j.at("smokeTextureFile").get<std::string>();
            s.plumeDT = j.at("plumeDT").get<float>();
            s.subspheres = j.at("subspheres").get<int>();
            s.subsubspheres = j.at("subsubspheres").get<int>();
            s.initialPlumeSpeed = j.at("initialPlumeSpeed").get<float>();
            s.initialPlumeDensity = j.at("initialPlumeDensity").get<float>();
            s.ventRay = j.at("ventRay").get<float>();
            s.ventAltitude = j.at("ventAltitude").get<float>();
            s.initialEjectionAngle = j.at("initialEjectionAngle").get<float>();
        }
    }
}
