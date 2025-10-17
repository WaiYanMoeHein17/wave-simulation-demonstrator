#ifndef DIFFERENCE_HPP
#define DIFFERENCE_HPP

#include <filesystem>
#include <iostream>

#include "Screen.hpp"
#include <opencv2/opencv.hpp>


namespace fs = std::filesystem;

class Difference : public Screen {
private:
    cv::Mat color;
    cv::Mat gray;
    int counter = 0;
public:
    Difference(std::string colorPath) : Screen() {
        color = cv::imread(colorPath, cv::IMREAD_COLOR);
        if (color.empty()) {
            std::cerr << "Error: Could not load image at " << colorPath << std::endl;
            return;
        }

        std::string grayPath = colorPath;
        size_t pos = grayPath.find("_col");
        if (pos != std::string::npos) {
            grayPath.replace(pos, 4, "_gray");
        }

        gray = cv::imread(grayPath, cv::IMREAD_GRAYSCALE);
        if (gray.empty()) {
            std::cerr << "Error: Could not load image at " << grayPath << std::endl;
        }

        update(color);
    }
    void toggle(const cv::Mat& terrainDepth) {
        if (color.size() != terrainDepth.size()) {
            cv::resize(color, color, terrainDepth.size(), 0, 0, cv::INTER_CUBIC);
        }
        if (gray.size() != terrainDepth.size()) {
            cv::resize(gray, gray, terrainDepth.size(), 0, 0, cv::INTER_CUBIC);
        }
        if(counter++ > 2) counter = 0;
        if (counter == 0) {
            update(color);
        } else if (counter == 1) {
            update(gray);
        } else if (counter == 2) {
            cv::Mat diff;
            cv::subtract(terrainDepth, gray, diff, cv::noArray(), CV_32F);
            cv::normalize(diff, diff, 0, 255, cv::NORM_MINMAX);
            diff.convertTo(diff, CV_8UC1);
            cv::Mat diffColor;
            cv::applyColorMap(diff, diffColor, cv::COLORMAP_VIRIDIS);

            update(diffColor);
        }
    }
};

#endif // DIFFERENCE_HPP