#ifndef WATER_HPP
#define WATER_HPP

#include <glad/glad.h>
#include <opencv2/opencv.hpp>

class Water {
private:
public:
    GLuint texture;
    cv::Mat waterDepthMat;
    cv::Mat waterTextureMat;
    float maxValue = 0.0f;
    Water(int width, int height) : texture(0) {
        waterDepthMat = cv::Mat::zeros(height, width, CV_32F);
        // cv::imshow("Water", waterDepthMat);
        // cv::setMouseCallback("Water", onMouse, this);
    };

    static void onMouse(int event, int x, int y, int flags, void* userdata)
    {
        if (event == cv::EVENT_LBUTTONDOWN || (event == cv::EVENT_MOUSEMOVE && (flags & cv::EVENT_FLAG_LBUTTON))) {
            Water* water = reinterpret_cast<Water*>(userdata);
            if (water) {
                int brushSize = 168;
                float sigma = 24.0f;
                float alpha = 0.1f;
                water->drawGaussianBrush(x, y, brushSize, sigma, alpha);
            }
        }
    }

    void drawGaussianBrush(int x, int y, int size, float sigma, float alpha)
    {
        int halfSize = size / 2;

        // Generate Gaussian kernel using OpenCV function
        cv::Mat kernelX = cv::getGaussianKernel(size, sigma, CV_32F);
        cv::Mat kernelY = cv::getGaussianKernel(size, sigma, CV_32F);
        cv::Mat gaussianKernel = kernelX * kernelY.t();

        // Normalize and scale by alpha
        double maxVal;
        cv::minMaxLoc(gaussianKernel, nullptr, &maxVal);
        gaussianKernel = (gaussianKernel / maxVal) * alpha;

        // apply the kernel to waterDepthMat
        for (int i = 0; i < size; ++i) {
            int posY = y + i - halfSize;
            if (posY < 0 || posY >= waterDepthMat.rows) continue;

            for (int j = 0; j < size; ++j) {
                int posX = x + j - halfSize;
                if (posX < 0 || posX >= waterDepthMat.cols) continue;

                waterDepthMat.at<float>(posY, posX) += gaussianKernel.at<float>(i, j);
                // waterDepthMat.at<float>(posY, posX) = std::min(2.55f * 2, std::max(0.0f, waterDepthMat.at<float>(posY, posX) + gaussianKernel.at<float>(i, j)));
            }
        }
        toGL();
        // cv::imshow("Water", waterDepthMat);
    }

    void clear()
    {
        waterDepthMat.setTo(0);
        toGL();
    }
    void increase()
    {
        waterDepthMat += 0.05f;
        toGL();
    }
    void decrease()
    {
        waterDepthMat -= 0.05f;
        toGL();
    }

    void toGL()
    {
        if (texture == 0) {
            glGenTextures(1, &texture);
        }
        glBindTexture(GL_TEXTURE_2D, texture);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        
        // Update maxValue to track the maximum value in waterDepthMat
        // double max;
        // cv::minMaxLoc(waterDepthMat, nullptr, &max);
        // maxValue = static_cast<float>(max);

        // CV_32F -> CV_8UC1
        cv::Mat normalizedMat;
        waterDepthMat.convertTo(normalizedMat, CV_8UC1, 255.0);

        // jet
        // cv::applyColorMap(normalizedMat, jetMat, cv::COLORMAP_JET);
        cv::applyColorMap(normalizedMat, waterTextureMat, cv::COLORMAP_OCEAN);

        // Clamp waterDepthMat values between 0.0 and 3.0
        cv::threshold(waterDepthMat, waterDepthMat, 1.0, 1.0, cv::THRESH_TRUNC);
        cv::threshold(waterDepthMat, waterDepthMat, 0.0, 0.0, cv::THRESH_TOZERO);
        
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, jetMat.cols, jetMat.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, jetMat.data);
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, normalizedMat.cols, normalizedMat.rows, 0, GL_RED, GL_UNSIGNED_BYTE, normalizedMat.data);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, waterDepthMat.cols, waterDepthMat.rows, 0, GL_RED, GL_FLOAT, waterDepthMat.data);
        


        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, rgbaMat.cols, rgbaMat.rows, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgbaMat.data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glBindTexture(GL_TEXTURE_2D, 0);
    }
    void saveDepth(const std::string& filename)
    {
        // Save as 32F EXR or TIFF for lossless float storage
        cv::imwrite(filename, waterDepthMat);
    }
    void loadDepth(const std::string& filename)
    {
        // Load as-is (should be 32F if saved as EXR/TIFF)
        waterDepthMat = cv::imread(filename, cv::IMREAD_UNCHANGED);
        if (waterDepthMat.empty()) {
            std::cerr << "Failed to load depth image: " << filename << std::endl;
            return;
        }
        if (waterDepthMat.type() != CV_32F) {
            std::cerr << "Failed to load depth image: " << filename << std::endl;
            return;
        }
        toGL();
    }
    
};

#endif // WATER_HPP