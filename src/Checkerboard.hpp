#ifndef CHECKERBOARD_HPP
#define CHECKERBOARD_HPP

#include "Screen.hpp"
#include <opencv2/opencv.hpp>
#include <vector>

class Checkerboard : public Screen {
public:
    std::vector<cv::Point> points;

    Checkerboard(int width, int height)
        : Screen()
    {
        cv::Mat checkerboard = createCheckerboard(width, height, points);
        update(checkerboard);
    }
    void debugCheckerboard() {
        for (const auto& point : this->points) {
            std::cout << "Point: (" << point.x << ", " << point.y << ")" << std::endl;
            cv::circle(image, point, 5, cv::Scalar(0, 0, 255), -1);
        }
        cv::imshow("Checkerboard", this->image);
    }
    void reset() {
        points.clear();
        cv::Mat checkerboard = createCheckerboard(image.cols, image.rows, points);
        update(checkerboard);
    }

private:
    static cv::Mat createCheckerboard(int width, int height, std::vector<cv::Point>& points) {
        cv::Mat checkerboardMat = cv::Mat::zeros(height, width, CV_8UC3);
        checkerboardMat.setTo(cv::Scalar(255, 255, 255));

        int squareSize = 80;
        int checkerboardSize = 8;
        int startX = (width - checkerboardSize * squareSize) / 2;
        int startY = (height - checkerboardSize * squareSize) / 2;

        for (int i = 0; i < checkerboardSize; ++i) {
            for (int j = 0; j < checkerboardSize; ++j) {
                cv::Rect square(startX + j * squareSize, startY + i * squareSize, squareSize, squareSize);
                cv::Scalar color = (i + j) % 2 == 0 ? cv::Scalar(0, 0, 0) : cv::Scalar(255, 255, 255);
                cv::rectangle(checkerboardMat, square, color, cv::FILLED);

                if (i < checkerboardSize - 1 && j < checkerboardSize - 1) {
                    points.emplace_back(startX + (j+1) * squareSize, 
                                        startY + (i+1) * squareSize);
                }
            }
        }

        return checkerboardMat;
    }

};

#endif // CHECKERBOARD_HPP
