#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <iostream>
#include <string>

#include <GLFW/glfw3.h>
#include <opencv2/opencv.hpp>
#include "Visualisation.hpp"
#include "Simulation.hpp"
#include "Water.hpp"
#include "ColorMap.hpp"
#include "Difference.hpp"

struct WindowState {
    Visualisation *visualisation;
    bool debugWindows;

    Simulation *simulation;
    
    Water *waterCanvas;

    Difference *differenceCanvas;

    bool isLeftDragging;
    bool isRightDragging;

    bool enableFilter = true;
    bool saveNext = false;
    bool shiftPressed = false;

    std::vector<cv::Point> markers;

    WindowState() : visualisation(nullptr), debugWindows(false), isLeftDragging(false), isRightDragging(false) {}
};

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    WindowState* state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));

    if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
    if (key == GLFW_KEY_T && action == GLFW_PRESS) {
        std::cout << "T key pressed" << std::endl;
    }

    if (key == GLFW_KEY_F && action == GLFW_PRESS) {
        state->visualisation->toggleGradientColor();
        std::cout << "Gradient Color: " << (state->visualisation->useGradientColor ? "ON" : "OFF") << std::endl;
    }
    if (key == GLFW_KEY_G && action == GLFW_PRESS) {
        state->visualisation->toggleGrayscale();
        std::cout << "Grayscale: " << (state->visualisation->useGrayscale ? "ON" : "OFF") << std::endl;
    }

    if (key == GLFW_KEY_RIGHT_BRACKET && action == GLFW_PRESS) {
        state->visualisation->incrementContourLineFactor();
        std::cout << "Contour Line Factor: " << state->visualisation->contourLineFactor << std::endl;
    }
    if (key == GLFW_KEY_LEFT_BRACKET && action == GLFW_PRESS) {
        state->visualisation->decrementContourLineFactor();
        std::cout << "Contour Line Factor: " << state->visualisation->contourLineFactor << std::endl;
    }
    if (key == GLFW_KEY_P && action == GLFW_PRESS) {
        state->visualisation->resetContourLineFactor();
        std::cout << "Contour Line Factor reset to: " << state->visualisation->contourLineFactor << std::endl;
    }
    if (key == GLFW_KEY_O && action == GLFW_PRESS) {
        state->visualisation->setContourLineFactor(0.0f);
        std::cout << "Contour Line Factor set to: " << state->visualisation->contourLineFactor << std::endl;
    }
    if (key == GLFW_KEY_GRAVE_ACCENT && action == GLFW_PRESS) {
        state->debugWindows = !state->debugWindows;
        std::cout << "Debug Windows: " << (state->debugWindows ? "ON" : "OFF") << std::endl;
        if (!state->debugWindows) {
            cv::destroyAllWindows();
        }
    }
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        state->visualisation->saveImage("image.png");
    }

    if (key == GLFW_KEY_E && action == GLFW_PRESS) {
        state->simulation->loadSequencePaths();
    }

    if (key == GLFW_KEY_C && action == GLFW_PRESS) {
        std::cout << "Save NetCDF" << std::endl;
        state->simulation->saveToNetCDF(state->visualisation->terrainImage);
    }

    if (key == GLFW_KEY_I && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
        state->visualisation->simulationOffset += 0.025f;
        std::cout << "Simulation Offset increased to: " << state->visualisation->simulationOffset << std::endl;
    }
    if (key == GLFW_KEY_U && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
        state->visualisation->simulationOffset -= 0.025f;
        std::cout << "Simulation Offset decreased to: " << state->visualisation->simulationOffset << std::endl;
    }
    // if (key == GLFW_KEY_K && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
    //     state->visualisation->simulationScale *= 1.025f;
    //     std::cout << "Simulation Scale increased to: " << state->visualisation->simulationScale << std::endl;
    // }
    // if (key == GLFW_KEY_J && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
    //     state->visualisation->simulationScale /= 1.025f;
    //     std::cout << "Simulation Scale decreased to: " << state->visualisation->simulationScale << std::endl;
    // }
    if (key == GLFW_KEY_T && action == GLFW_PRESS) {
        state->visualisation->toggleWaterTexture();
        std::cout << "Water Texture: " << (state->visualisation->useWaterTexture ? "ON" : "OFF") << std::endl;
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS) {
        state->visualisation->togglePause();
        std::cout << "Visualisation Pause: " << (state->visualisation->isPaused() ? "ON" : "OFF") << std::endl;
    }
    if (key == GLFW_KEY_SLASH && action == GLFW_PRESS) {
        std::cout << "--------------------------" << std::endl;
        std::cout << "Running Simulation" << std::endl;
        // if (state->debugWindows) {
        //     cv::imshow("Water Canvas", state->waterCanvas->waterDepthMat);
        //     cv::imshow("Terrain Image", state->visualisation->terrainImage);
        //     cv::waitKey(0);
        // }
        state->visualisation->paused = true;
        state->simulation->saveToNetCDF(state->visualisation->terrainImage);
        state->simulation->saveToNetCDF(state->waterCanvas->waterDepthMat, "water.nc");
        state->simulation->testLoadSequence();
        state->visualisation->simulationOffset = state->simulation->pastOffset;
    }
    if (key == GLFW_KEY_W && action == GLFW_PRESS) {
        state->simulation->reset();
        state->visualisation->simulationOffset = 0.0f;
        std::cout << "Simulation reset" << std::endl;
    }
    if (key == GLFW_KEY_H && action == GLFW_PRESS) {
        state->visualisation->flipY();
        std::cout << "Flipped Y" << std::endl;
    }
    if (key == GLFW_KEY_V && action == GLFW_PRESS) {
        state->visualisation->flipX();
        std::cout << "Flipped X" << std::endl;
    }
    if (action == GLFW_PRESS && (key == GLFW_KEY_1 || key == GLFW_KEY_2 || key == GLFW_KEY_3 || key == GLFW_KEY_4)) {
        state->visualisation->updateVertexFromMouseInput(window, key, action);
    }
    // if(key == GLFW_KEY_W && action == GLFW_PRESS) {
    //     cv::imshow("Water Canvas", state->waterCanvas->waterMat);
    //     cv::waitKey(1);
    //     state->simulation->saveToNetCDF(state->waterCanvas->waterMat, "water.nc");
    // }


    // Water Canvas Controls
    if (key == GLFW_KEY_EQUAL && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
        state->waterCanvas->increase();
        std::cout << "Water Canvas increased" << std::endl;
    }
    if (key == GLFW_KEY_MINUS && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
        state->waterCanvas->decrease();
        std::cout << "Water Canvas decreased" << std::endl;
    }
    if (key == GLFW_KEY_0 && action == GLFW_PRESS) {
        state->waterCanvas->clear();
        std::cout << "Water Canvas cleared" << std::endl;
    }
    if (key == GLFW_KEY_SEMICOLON && action == GLFW_PRESS) {
        state->visualisation->colorMapIndex = (state->visualisation->colorMapIndex + 1) % COLOR_MAPS.size();
        std::cout << "Colormap Index increased to: " << state->visualisation->colorMapIndex << std::endl;
    }
    if (key == GLFW_KEY_D && action == GLFW_PRESS) {
        state->simulation->toggleField();
        std::cout << "Simulation Field: " << (state->simulation->isFieldEnabled() ? "ON" : "OFF") << std::endl;
    }
    if(key == GLFW_KEY_F6 && action == GLFW_PRESS) {
        // Save water depth with a timestamp in ../waters
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << "../waters/waterDepth_"
            << std::put_time(&tm, "%Y%m%d_%H%M%S")
            << ".tiff";
        state->waterCanvas->saveDepth(oss.str());
        std::cout << "Water Canvas saved" << std::endl;
    }
    if(key == GLFW_KEY_PERIOD && action == GLFW_PRESS) {
        // state->waterCanvas->loadDepth("waterDepth.tiff");
        static int tiffIndex = 0;
        std::vector<std::string> tiffFiles;
        if (tiffFiles.empty()) {
            std::string dir = "../waters";
            cv::glob(dir + "/*.tiff", tiffFiles, false);
            std::sort(tiffFiles.begin(), tiffFiles.end());
        }
        if (!tiffFiles.empty()) {
            tiffIndex = (tiffIndex + 1) % tiffFiles.size();
            state->waterCanvas->loadDepth(tiffFiles[tiffIndex]);
            std::cout << "Loaded: " << tiffFiles[tiffIndex] << std::endl;
        } else {
            std::cout << "No .tiff files found in ../waters" << std::endl;
        }
    }

    if (key == GLFW_KEY_F1 && action == GLFW_PRESS) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        cv::Point newPoint(static_cast<int>(xpos), static_cast<int>(ypos));
        if (state->markers.size() >= 2) {
            state->markers.erase(state->markers.begin());
        }
        state->markers.push_back(newPoint);
        std::cout << "Added point: (" << newPoint.x << ", " << newPoint.y << ")" << std::endl;
    }
    if (key == GLFW_KEY_F2 && action == GLFW_PRESS) {
        if (state->markers.size() == 2) {
            int windowWidth, windowHeight;
            glfwGetWindowSize(window, &windowWidth, &windowHeight);
            state->visualisation->samplePointsAlongLine(state->markers[0], state->markers[1], 100, 
                windowWidth, windowHeight,
            state->waterCanvas->waterDepthMat);
        }
    }
    if (key == GLFW_KEY_F3 && action == GLFW_PRESS) {
        cv::Mat glImage = state->visualisation->getGLImage();
        cv::imshow("GL Image", glImage);
        cv::waitKey(1);
    }

    if (key == GLFW_KEY_F4 && action == GLFW_PRESS) {
        state->enableFilter = !state->enableFilter;
        std::cout << "Filter: " << (state->enableFilter ? "ON" : "OFF") << std::endl;
    }
    if (key == GLFW_KEY_F5 && action == GLFW_PRESS) {
        state->saveNext = true;
        std::cout << "Save next frame" << std::endl;
    }
    if (key == GLFW_KEY_Z && action == GLFW_PRESS) {
        state->differenceCanvas->toggle(state->visualisation->terrainImage);
    }
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)     state->shiftPressed = true;
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)   state->shiftPressed = false;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    WindowState* state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        state->isLeftDragging = true;
        std::cout << "Left mouse button pressed at position (" << xpos << ", " << ypos << ")" << std::endl;
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        state->isLeftDragging = false;
        std::cout << "Left mouse button released at position (" << xpos << ", " << ypos << ")" << std::endl;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        state->isRightDragging = true;
        std::cout << "Right mouse button pressed at position (" << xpos << ", " << ypos << ")" << std::endl;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
        state->isRightDragging = false;
        std::cout << "Right mouse button released at position (" << xpos << ", " << ypos << ")" << std::endl;
    }
}

void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
    WindowState* state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
    if (state->isLeftDragging) {
        float strength = state->shiftPressed ? 0.1f : 0.05f;
        // std::cout << "Left mouse dragging at position (" << xpos << ", " << ypos << ")" << std::endl;
        state->waterCanvas->drawGaussianBrush(static_cast<int>(xpos), static_cast<int>(ypos), 168, 24.0f, strength);
    }
    if (state->isRightDragging) {
        float strength = state->shiftPressed ? 0.1f : 0.05f;
        // std::cout << "Right mouse dragging at position (" << xpos << ", " << ypos << ")" << std::endl;
        state->waterCanvas->drawGaussianBrush(static_cast<int>(xpos), static_cast<int>(ypos), 168, 24.0f, -strength);
    }
}

#endif // WINDOW_HPP