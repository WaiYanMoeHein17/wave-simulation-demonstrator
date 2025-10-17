#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include <vector>
#include <numeric>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>

#include <opencv2/opencv.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <umbridge/umbridge.h>
#include <netcdf>

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>

#define SIMULATION_WIDTH 800
#define SIMULATION_HEIGHT 600

namespace fs = std::filesystem;

class Simulation {
private:
fs::path inputPath;
fs::path outputPath;
std::string host;

std::thread simulationThread;
std::mutex sequenceMutex;
std::vector<fs::path> sequencePaths;
std::vector<cv::Mat> frames;
std::vector<cv::Mat> textureFrames;
std::vector<float> maxDepthValues;

bool first = true;
std::vector<std::pair<int, int>> normalizedCoordinates;
bool field = false;

public:
    double simStartTime = 0.0;
    double pastOffset = 0.0;
    unsigned int currentFrame = 0;
    GLuint texture = 0;
    std::atomic<bool> isRunning;

    Simulation(std::string inPath, std::string outPath, std::string hostAddress = "http://localhost:4242") 
        : inputPath(inPath), outputPath(outPath), host(hostAddress),
        isRunning(false) {}

    void toggleField() {
        field = !field;
    }
    bool isFieldEnabled() {
        return field;
    }

    cv::Mat *getCurrentTexture() {
        std::lock_guard<std::mutex> lock(sequenceMutex);
        if (currentFrame - 1 >= sequencePaths.size()) {
            return nullptr;
        }
        return &textureFrames[currentFrame - 1];
    }

    cv::Mat *getCurrentFrame() {
        std::lock_guard<std::mutex> lock(sequenceMutex);
        if (currentFrame - 1 >= sequencePaths.size()) {
            return nullptr;
        }
        return &frames[currentFrame - 1];
    }

    std::vector<fs::path> getSequencePaths(const std::string &path) {
        fs::path dir(path);
        if (!fs::exists(dir) || !fs::is_directory(dir)) {
            std::cerr << "Invalid directory path: " << path << std::endl;
        }
    
        std::vector<fs::path> files;
        for (const auto &entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".vtk") {
                files.push_back(entry.path());
            }
        }
    
        std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
            std::string numA = a.stem().string().substr(12, a.stem().string().length() - 12 - 7);
            std::string numB = b.stem().string().substr(12, b.stem().string().length() - 12 - 7);
            return std::stoi(numA) < std::stoi(numB);
        });
    
        return files;
    }

    void reset() {
        currentFrame = 0;
        frames.clear();
        normalizedCoordinates.clear();
        sequencePaths.clear();
        maxDepthValues.clear();
        textureFrames.clear();
        first = true;
    }

    void triggerSimulation() {}
    void testLoadSequence() {
        if (simulationThread.joinable()) {
            std::cerr << "Simulation thread is already running. Please wait for it to complete." << std::endl;
            return;
        }
        // reset
        reset();

        simulationThread = std::thread([this]() {
            isRunning = true;
            try {
                std::cout << "Triggering simulation [" << host << "]" << std::endl;
                umbridge::HTTPModel client(host, "forward");
                std::vector<std::vector<double>> inputs {{100.0, 100.0}};
                std::vector<std::vector<double>> outputs = client.Evaluate(inputs);
                std::cout << "--------------------------------" << std::endl;
                std::cout << "Simulation Output: ";
                for (const auto& value : outputs[0]) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;

                std::lock_guard<std::mutex> lock(sequenceMutex);
                sequencePaths = getSequencePaths(inputPath.string());
                std::cout << "Found " << sequencePaths.size() << " sequence files." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "An error occurred during simulation: " << e.what() << std::endl;
            }
            isRunning = false;
        });
        simulationThread.detach();
    }

    void loadSequencePaths() {
        reset();
        std::lock_guard<std::mutex> lock(sequenceMutex);
        sequencePaths = getSequencePaths(inputPath.string());
        std::cout << "Found " << sequencePaths.size() << " sequence files." << std::endl;
    }

    unsigned int frameCount() {
        std::lock_guard<std::mutex> lock(sequenceMutex);
        return sequencePaths.size();
    }

    void advanceFrame(float &maxValue) {
        std::lock_guard<std::mutex> lock(sequenceMutex);
        if (sequencePaths.empty()) {
            // std::cerr << "No sequence paths available." << std::endl;
            return;
        }
        fs::path nextPath = sequencePaths[currentFrame++];
        if (currentFrame >= sequencePaths.size()) {
            currentFrame = 0;
        }
        // std::cout << "Loading frame: " << nextPath.stem().string() << std::endl;

        if (frames.size() == sequencePaths.size()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(80));
            maxValue = maxDepthValues[currentFrame];
            toGL(frames[currentFrame]);
            return;
        }


        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->ReadAllScalarsOn();
        reader->SetFileName(nextPath.string().c_str());
        reader->Update();

        vtkSmartPointer<vtkUnstructuredGrid> ugrid = reader->GetOutput();
        if (!ugrid) {
            std::cerr << "Failed to read the VTK file." << std::endl;
            return;
        }

        vtkDataArray* scalarQ = ugrid->GetCellData()->GetArray("Q");
        if (!scalarQ) {
            std::cerr << "Scalar array 'Q' not found!" << std::endl;
            return;
        }

    
        // int numTuples = scalarQ->GetNumberOfTuples();
        int numCells = ugrid->GetNumberOfCells();
        // int numComponents = scalarQ->GetNumberOfComponents();

        if(first) {
            double minX = std::numeric_limits<double>::max();
            double minY = std::numeric_limits<double>::max();
            double maxX = std::numeric_limits<double>::lowest();
            double maxY = std::numeric_limits<double>::lowest();

            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

            for (int i = 0; i < numCells; ++i) {
                vtkCell* cell = ugrid->GetCell(i);
                if (!cell) continue;

                vtkSmartPointer<vtkPoints> cellPoints = cell->GetPoints();
                if (!cellPoints) continue;

                // cell center
                double centerX = 0, centerY = 0, centerZ = 0;
                int numCellPoints = cellPoints->GetNumberOfPoints();
                for (int j = 0; j < numCellPoints; ++j) {
                    double p[3];
                    cellPoints->GetPoint(j, p);
                    centerX += p[0];
                    centerY += p[1];
                    centerZ += p[2];
                }
                centerX /= numCellPoints;
                centerY /= numCellPoints;
                centerZ /= numCellPoints;

                // update bounding box
                minX = std::min(minX, centerX);
                minY = std::min(minY, centerY);
                maxX = std::max(maxX, centerX);
                maxY = std::max(maxY, centerY);

                points->InsertNextPoint(centerX, centerY, centerZ);
            }
            // norm
            for (int i = 0; i < points->GetNumberOfPoints(); ++i) {
                double p[3];
                points->GetPoint(i, p);
                int x = static_cast<int>((p[0] - minX) / (maxX - minX) * SIMULATION_WIDTH - 1);
                int y = static_cast<int>((p[1] - minY) / (maxY - minY) * SIMULATION_HEIGHT - 1);
                normalizedCoordinates.emplace_back(x, y);
            }

            // std::cout << "Bounding box: (" << minX << ", " << minY << ") to (" << maxX << ", " << maxY << ")" << std::endl;
            first = false;
        }

        // cv::Mat depthMap = cv::Mat::zeros(SIMULATION_HEIGHT, SIMULATION_WIDTH, CV_32FC1);
        // for (size_t i = 0; i < normalizedCoordinates.size(); ++i) {
        //     const auto& [x, y] = normalizedCoordinates[i];
        //     depthMap.at<float>(y, x) = static_cast<float>(scalarQ->GetComponent(i, 0));
        // }

        //  .at<>() ===== error checking
        cv::Mat depthMap = cv::Mat::zeros(SIMULATION_HEIGHT, SIMULATION_WIDTH, CV_32FC1);
        for (size_t i = 0; i < normalizedCoordinates.size(); ++i) {
            // const auto& [x, y] = normalizedCoordinates[i];
            // depthMap.ptr<float>(y)[x] = static_cast<float>(scalarQ->GetComponent(i, 0));
            const auto& [x, y] = normalizedCoordinates[i];
            if (x >= 0 && x < SIMULATION_WIDTH && y >= 0 && y < SIMULATION_HEIGHT) {
                depthMap.ptr<float>(y)[x] = static_cast<float>(scalarQ->GetComponent(i, 0));
            }
        }

        double minDepth, maxDepth;
        cv::minMaxLoc(depthMap, &minDepth, &maxDepth);
        std::cout << "Depth Map - Min: " << minDepth << ", Max: " << maxDepth << std::endl;
        maxDepthValues.push_back(maxDepth);
        maxValue = maxDepth;

        // heavy gaussian for complete image
        cv::Mat interpolatedMap;
        cv::GaussianBlur(depthMap, interpolatedMap, cv::Size(61, 61), 0);
        cv::normalize(interpolatedMap, interpolatedMap, 0.0f, 1.0f, cv::NORM_MINMAX);
        // cv::imshow("INMAP", interpolatedMap);
        // cv::waitKey(1);
        // cv::GaussianBlur(depthMap, interpolatedMap, cv::Size(5, 5), 0);

        // norm
        cv::Mat colorMap;
        cv::Mat normalizedMap = interpolatedMap.clone();
        cv::normalize(normalizedMap, normalizedMap, 0, 255, cv::NORM_MINMAX);
        normalizedMap.convertTo(normalizedMap, CV_8U);
        cv::applyColorMap(normalizedMap, colorMap, cv::COLORMAP_JET);

        if (field) {
            colorMap.setTo(cv::Scalar(255, 255, 255));
        }
        
        for (size_t i = 0; i < normalizedCoordinates.size(); i += 20) {
            const auto& [x, y] = normalizedCoordinates[i];
            float xMag = static_cast<float>(scalarQ->GetComponent(i, 1));
            float yMag = static_cast<float>(scalarQ->GetComponent(i, 2));

            float direction = atan2(yMag, xMag);
            float magnitude = sqrt(xMag * xMag + yMag * yMag) * 200;
            if (magnitude < 1.0) continue;
            cv::Point arrowStart(x, y);
            cv::Point arrowEnd(x + static_cast<int>(magnitude * cos(direction)), y + static_cast<int>(magnitude * sin(direction)));
            // cv::arrowedLine(colorMap, arrowStart, arrowEnd, field ? cv::Scalar(0, 0, 0) : cv::Scalar(255, 255, 255), 1, cv::LINE_AA, 0, 0.5);
            cv::arrowedLine(colorMap, arrowStart, arrowEnd, field ? cv::Scalar(0, 0, 0) : cv::Scalar(255, 255, 255), 2, cv::LINE_AA, 0, 0.5);
        }
        textureFrames.push_back(colorMap);
        // cv::normalize(interpolatedMap, interpolatedMap, 0, 255, cv::NORM_MINMAX);
        // interpolatedMap.convertTo(interpolatedMap, CV_8U);

        // cv::Mat colorMap;
        // cv::applyColorMap(interpolatedMap, colorMap, cv::COLORMAP_JET);
        
        
        // for (size_t i = 0; i < normalizedCoordinates.size(); i += 20) {
        //     const auto& [x, y] = normalizedCoordinates[i];
        //     float xMag = static_cast<float>(scalarQ->GetComponent(i, 1));
        //     float yMag = static_cast<float>(scalarQ->GetComponent(i, 2));

        //     float direction = atan2(yMag, xMag);
        //     float magnitude = sqrt(xMag * xMag + yMag * yMag) * 100;
        //     if (magnitude < 1.0) continue;
        //     cv::Point arrowStart(x, y);
        //     cv::Point arrowEnd(x + static_cast<int>(magnitude * cos(direction)), y + static_cast<int>(magnitude * sin(direction)));
        //     cv::arrowedLine(colorMap, arrowStart, arrowEnd, cv::Scalar(255, 255, 255), 1, cv::LINE_AA, 0, 0.5);
        // }

        frames.push_back(interpolatedMap);
        toGL(interpolatedMap);

        // frames.push_back(colorMap);
        // toGL(colorMap);
        // cv::imshow("Depth Map", colorMap);
        // cv::waitKey(1);
    }
    void toGL(const cv::Mat& frame)
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
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, flipped.cols, flipped.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, flipped.data);

        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frame.cols, frame.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, frame.data);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, frame.cols, frame.rows, 0, GL_RED, GL_FLOAT, frame.data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void saveToNetCDF(cv::Mat image, std::string filename="output.nc")
    {
        cv::resize(image, image, cv::Size(800, 600), 0, 0, cv::INTER_LINEAR);
        size_t ny = image.rows;
        size_t nx = image.cols;

        double minVal, maxVal;
        cv::minMaxLoc(image, &minVal, &maxVal);
        std::cout << "min: " << minVal << ", max: " << maxVal << std::endl;

        if (filename == "water.nc") {
            this->pastOffset = minVal * 2.55;
            std::cout << "pastOffset: " << this->pastOffset << std::endl;
        }

        std::vector<float> x(nx), y(ny);
        std::iota(x.begin(), x.end(), 0.0f);
        std::iota(y.begin(), y.end(), 0.0f);

        fs::path heightmap_path = fs::path(outputPath) / filename;
        netCDF::NcFile ncFile(heightmap_path.string(), netCDF::NcFile::replace);

        netCDF::NcDim xDim = ncFile.addDim("x", nx);
        netCDF::NcDim yDim = ncFile.addDim("y", ny);

        netCDF::NcVar xVar = ncFile.addVar("x", netCDF::ncFloat, xDim);
        netCDF::NcVar yVar = ncFile.addVar("y", netCDF::ncFloat, yDim);
        netCDF::NcVar zVar = ncFile.addVar("elevation", netCDF::ncFloat, {yDim, xDim});

        std::vector<float> zData(ny * nx);
        // const uchar* imageData = image.ptr<uchar>();
        // std::transform(imageData, imageData + (ny * nx), zData.begin(), [](uchar val) {
        //     return static_cast<float>(val);
        // });

        // !different image types
        if (image.type() == CV_8U) {
            const uchar* imageData = image.ptr<uchar>();
            std::transform(imageData, imageData + (ny * nx), zData.begin(), [](uchar val) {
                return static_cast<float>(val);
            });
        } else if (image.type() == CV_32F) {
            const float* imageData = image.ptr<float>();
            std::transform(imageData, imageData + (ny * nx), zData.begin(), [](float val) {
                return static_cast<float>(val) * 255;
            });
        } else {
            std::cerr << "Unsupported image type!" << std::endl;
            return;
        }
        
        
        // std::vector<float> zData(ny * nx);
        // for (size_t i = 0; i < ny; ++i) {
        //     for (size_t j = 0; j < nx; ++j) {
        //         zData[i * nx + j] = static_cast<float>(image.at<uchar>(i, j));
        //     }
        // }

        xVar.putVar(x.data());
        yVar.putVar(y.data());
        zVar.putVar(zData.data());

        zVar.putAtt("units", "meters");
        zVar.putAtt("standard_name", "surface_elevation");

        std::cout << filename <<" NetCDF file created successfully!" << std::endl;
    }
};

#endif // SIMULATION_HPP
