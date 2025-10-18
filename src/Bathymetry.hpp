#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

#include <vector> 
#include <string> 
#include <memory> 
#include <iostream> 
#include <glm/glm.hpp>

struct CADModel {
    std::string name; 
    std::string filePath;
    glm::vec3 position;
    glm::vec3 scale;
    glm::vec3 rotation;
    float depth; 
    bool isVisible; 
    bool isSubmerged;

    std::vector<glm::evec3> vertices;
    std::vector<unsigned int> indices;
    std::vector<glm::uvec3> faces;

    CADModel() : position(0.0f), scale(1.0f), rotation(0.0f), depth(0.0f), isVisible(true), isSubmerged(false) {}
};

class Bathymetry {
    public: 
        Bathymetry(size_t gridWidth, size_t gridHeight); 
        ~Bathymetry();

        int addCADModel(const std::string& name, const std::string& filePath, const glm::vec3& position, const glm::vec3& scale, const glm::vec3& rotation, float depth);
        bool removeCADModel(int ModelID);
        void updateModelTransform(int modelID, const glm::vec3& position, const glm::vec3& scale, const glm::vec3& rotation);
        void setModelDepth(int modelID, float depth);

        // Load CAD model from file
        bool loadSTL(const std::string& filePath, CADModel& model);
        bool loadOBJ(const std::string& filePath, CADModel& model);
        bool loadPLY(const std::string& filePath, CADModel& model);

        // Combine camera depth with CAD bathymetry data
        void combineBathymetry(const float* cameraDepth, float* outputDepth); 

        // Rasterize CAD models to depth grid
        void rasterizeModels(); 
    
        const float* getBathymetryGrid() const {return bathymetryGrid.data(); }
        float* getBathymetryGrid() {return bathymetryGrid.data(); }

        // Visualization 
        void renderModels(); 
        void setGridResolution(size_t width, size_t height);

        // UI Helpers 
        const std::vector<CADModel>& getModels() const { return models; }
        void clearModels(); 

        bool saveConfig(const std::string& filePath);
        bool loadConfig(const std::string& filePath);

    private:
        size_t gridWidth;
        size_t gridHeight;
        alignas(32) std::vector<float> bathymetryGrid; 
        std::vector<CADModel> models; 

        void rasterizeTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, float depth);
        bool pointInTriangle(const glm::vec2& p, const glm::vec2& v0, const glm::vec2& v1, const glm::vec2& v2);
        float interpolateDepth(const glm::vec2& p, const glm::vec2& v0, const glm::vec2& v1, const glm::vec2& v2, float d0, float d1, float d2);

        glm::mat4 getModelMatrix(const CADModel& model);
        glm::vec3 worldToGrid(const glm::vec3& worldPos);
}; 

#endif // BATHYMETRY_HPP