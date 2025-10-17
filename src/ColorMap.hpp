#ifndef COLORMAP_HPP
#define COLORMAP_HPP

#include <vector>
#include <glm/glm.hpp>

std::vector<std::vector<glm::vec3>> COLOR_MAPS = {
    {
        glm::vec3(0.0f, 0.0f, 1.0f),  // Deep blue
        glm::vec3(0.0f, 0.5f, 1.0f),  // Light blue
        glm::vec3(0.0f, 1.0f, 0.0f),  // Green
        glm::vec3(1.0f, 1.0f, 0.0f),  // Yellow
        glm::vec3(1.0f, 0.5f, 0.0f),  // Orange
        glm::vec3(1.0f, 0.0f, 0.0f),  // Red
        glm::vec3(0.5f, 0.0f, 0.0f),  // Dark red
        glm::vec3(0.5f, 0.5f, 0.5f)   // Gray
    },
    {
        glm::vec3(0.0f, 0.0f, 1.0f),  // Deep blue
        glm::vec3(0.5f, 0.75f, 1.0f), // Light blue
        glm::vec3(0.0f, 1.0f, 0.0f),  // Green
        glm::vec3(0.5f, 1.0f, 0.5f),  // Light green
        glm::vec3(1.0f, 1.0f, 0.0f),  // Yellow
        glm::vec3(1.0f, 0.5f, 0.0f),  // Orange
        glm::vec3(1.0f, 0.0f, 0.0f),  // Red
        glm::vec3(0.5f, 0.0f, 0.0f)   // Dark red
    },
    {
        glm::vec3(0.2f, 0.0f, 0.0f),  // Dark red
        glm::vec3(0.6f, 0.0f, 0.0f),  // Red
        glm::vec3(1.0f, 0.3f, 0.0f),  // Orange-red
        glm::vec3(1.0f, 0.6f, 0.0f),  // Orange
        glm::vec3(1.0f, 0.8f, 0.2f),  // Yellow-orange
        glm::vec3(1.0f, 1.0f, 0.5f),  // Pale yellow
        glm::vec3(1.0f, 1.0f, 1.0f),  // White
        glm::vec3(0.5f, 0.5f, 0.5f)   // Gray
    },
    {
        glm::vec3(0.5f, 0.0f, 0.5f),  // Purple
        glm::vec3(0.75f, 0.0f, 1.0f), // Violet
        glm::vec3(1.0f, 0.0f, 1.0f),  // Magenta
        glm::vec3(1.0f, 0.5f, 1.0f),  // Light magenta
        glm::vec3(0.0f, 1.0f, 1.0f),  // Cyan
        glm::vec3(0.0f, 0.75f, 1.0f), // Light cyan
        glm::vec3(0.0f, 0.5f, 1.0f),  // Blue-cyan
        glm::vec3(0.0f, 0.25f, 0.5f)  // Deep blue-cyan
    },
    {
        glm::vec3(0.8f, 0.2f, 0.0f),  // Deep orange
        glm::vec3(1.0f, 0.5f, 0.0f),  // Orange
        glm::vec3(1.0f, 0.8f, 0.2f),  // Yellow-orange
        glm::vec3(1.0f, 1.0f, 0.5f),  // Pale yellow
        glm::vec3(0.8f, 0.6f, 1.0f),  // Lavender
        glm::vec3(0.5f, 0.3f, 1.0f),  // Purple
        glm::vec3(0.3f, 0.0f, 0.5f),  // Deep purple
        glm::vec3(0.2f, 0.0f, 0.3f)   // Dark violet
    },
    {
        glm::vec3(0.3f, 0.2f, 0.1f),  // Dark brown
        glm::vec3(0.5f, 0.3f, 0.2f),  // Brown
        glm::vec3(0.7f, 0.5f, 0.3f),  // Light brown
        glm::vec3(0.5f, 0.6f, 0.2f),  // Olive
        glm::vec3(0.2f, 0.8f, 0.2f),  // Green
        glm::vec3(0.0f, 0.6f, 0.0f),  // Dark green
        glm::vec3(0.0f, 0.4f, 0.0f),  // Deep green
        glm::vec3(0.1f, 0.3f, 0.1f)   // Forest green
    },
    {
        glm::vec3(0.3f, 0.2f, 0.1f),  // Dark brown
        glm::vec3(0.0f, 0.0f, 0.5f),  // Deep blue
        glm::vec3(0.5f, 0.75f, 1.0f), // Light blue
        glm::vec3(0.0f, 1.0f, 0.0f),  // Green
        glm::vec3(0.0f, 0.5f, 0.0f),  // Dark green
        glm::vec3(0.5f, 0.5f, 0.5f),  // Gray
        glm::vec3(1.0f, 1.0f, 1.0f),  // White
        glm::vec3(0.8f, 0.8f, 0.8f)   // Light gray
    },
    // ----------------------
    {
        glm::vec3(0.1f, 0.0f, 0.0f),  // Dark lava red
        glm::vec3(0.4f, 0.0f, 0.0f),  // Lava red
        glm::vec3(0.8f, 0.2f, 0.0f),  // Molten orange
        glm::vec3(1.0f, 0.5f, 0.0f),  // Bright lava orange
        glm::vec3(1.0f, 0.8f, 0.2f),  // Glowing yellow-orange
        glm::vec3(1.0f, 1.0f, 0.5f),  // Pale yellow
        glm::vec3(0.5f, 0.5f, 0.5f),  // Ash gray
        glm::vec3(0.2f, 0.2f, 0.2f)   // Dark ash
    }
};

#endif // COLORMAP_HPP