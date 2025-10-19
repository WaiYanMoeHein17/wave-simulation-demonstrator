#ifndef BATHYMETRY_EDITOR_HPP
#define BATHYMETRY_EDITOR_HPP

#include "Bathymetry.hpp"
#include <imgui.h>

class BathymetryEditor {
public:
    BathymetryEditor(Bathymetry& bathymetry) : bathymetry(bathymetry) {}
    
    void renderUI() {
        ImGui::Begin("Bathymetry Editor");
        
        if (ImGui::Button("Add Model")) {
            ImGui::OpenPopup("FileDialog");
        }
        
        // List existing models
        const auto& models = bathymetry.getModels();
        for (size_t i = 0; i < models.size(); ++i) {
            ImGui::PushID(i);
            
            if (ImGui::TreeNode(models[i].filepath.c_str())) {
                glm::vec3 pos = models[i].position;
                if (ImGui::DragFloat3("Position", &pos.x, 0.01f)) {
                    bathymetry.updateModelTransform(i, pos, models[i].rotation, models[i].scale);
                }
                
                glm::vec3 rot = models[i].rotation;
                if (ImGui::DragFloat3("Rotation", &rot.x, 1.0f)) {
                    bathymetry.updateModelTransform(i, models[i].position, rot, models[i].scale);
                }
                
                glm::vec3 scale = models[i].scale;
                if (ImGui::DragFloat3("Scale", &scale.x, 0.01f)) {
                    bathymetry.updateModelTransform(i, models[i].position, models[i].rotation, scale);
                }
                
                float depth = models[i].depth;
                if (ImGui::SliderFloat("Depth", &depth, -1.0f, 0.0f)) {
                    bathymetry.setModelDepth(i, depth);
                }
                
                if (ImGui::Button("Remove")) {
                    bathymetry.removeCADModel(i);
                }
                
                ImGui::TreePop();
            }
            
            ImGui::PopID();
        }
        
        ImGui::End();
    }
    
private:
    Bathymetry& bathymetry;
};

#endif