#ifndef REMOTE_HPP
#define REMOTE_HPP

#include <iostream>
#include <mutex>
#include <thread>
#include <fstream>

#include <httplib.h>

#include "Window.hpp"

class Remote {
    private:
    WindowState *state;
    public:
    std::mutex serverMutex;
    std::thread serverThread;
    void startServer() {
        httplib::Server svr;

        svr.Get("/", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received GET request on /" << std::endl;
            std::ifstream file("../remote/remote.html");
            if (file) {
                std::stringstream buffer;
                buffer << file.rdbuf();
                res.set_content(buffer.str(), "text/html");
            } else {
                res.status = 404;
                res.set_content("File not found", "text/plain");
            }
        });
        svr.Post("/toggle-greyscale", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /toggle-greyscale" << std::endl;
            state->visualisation->toggleGrayscale();
            res.set_content("Greyscale toggled", "text/plain");
        });
        svr.Post("/toggle-gradient", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /toggle-gradient" << std::endl;
            state->visualisation->toggleGradientColor();
            res.set_content("Gradient toggled", "text/plain");
        });
        svr.Post("/toggle-pause", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /toggle-pause" << std::endl;
            state->visualisation->togglePause();
            res.set_content("Pause toggled", "text/plain");
        });
        svr.Post("/increment-contour", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /increment-contour" << std::endl;
            state->visualisation->incrementContourLineFactor();
            res.set_content("Contour incremented", "text/plain");
        });
        svr.Post("/decrement-contour", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /decrement-contour" << std::endl;
            state->visualisation->decrementContourLineFactor();
            res.set_content("Contour decremented", "text/plain");
        });
        svr.Post("/reset-contour", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /reset-contour" << std::endl;
            state->visualisation->resetContourLineFactor();
            res.set_content("Contour reset", "text/plain");
        });
        svr.Post("/zero-contour", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /zero-contour" << std::endl;
            state->visualisation->setContourLineFactor(0.0f);
            res.set_content("Contour set to zero", "text/plain");
        });

        svr.Post("/save-image", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /save-image" << std::endl;
            state->visualisation->saveImage("image.png");
            res.set_content("Image saved", "text/plain");
        });

        svr.Post("/toggle-water", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /toggle-water" << std::endl;
            state->visualisation->toggleWaterTexture();
            res.set_content("Water texture toggled", "text/plain");
        });
        svr.Post("/toggle-field", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /toggle-field" << std::endl;
            state->simulation->toggleField();
            res.set_content("Field toggled", "text/plain");
        });

        svr.Post("/load-sequence", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /load-sequence" << std::endl;
            state->simulation->loadSequencePaths();
            res.set_content("Sequence loaded", "text/plain");
        });

        svr.Post("/increase-water-canvas", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /increase-water-canvas" << std::endl;
            state->waterCanvas->increase();
            res.set_content("Water Canvas increased", "text/plain");
        });

        svr.Post("/decrease-water-canvas", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /decrease-water-canvas" << std::endl;
            state->waterCanvas->decrease();
            res.set_content("Water Canvas decreased", "text/plain");
        });

        svr.Post("/clear-water-canvas", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /clear-water-canvas" << std::endl;
            state->waterCanvas->clear();
            res.set_content("Water Canvas cleared", "text/plain");
        });

        svr.Post("/next-colormap", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /next-colormap" << std::endl;
            state->visualisation->colorMapIndex = (state->visualisation->colorMapIndex + 1) % COLOR_MAPS.size();
            res.set_content("Colormap Index increased", "text/plain");
        });

        svr.Get("/terrain", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received GET request on /terrain" << std::endl;
            std::vector<uchar> buffer;
            cv::imencode(".png", state->visualisation->terrainImage, buffer);
            res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
        });
        svr.Get("/water", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received GET request on /water" << std::endl;
            std::vector<uchar> buffer;
            cv::Mat normalizedImage = state->waterCanvas->waterDepthMat.clone();
            normalizedImage.convertTo(normalizedImage, CV_8UC1, 255.0);
            // cv::normalize(normalizedImage, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_8UC1);
            
            cv::imencode(".png", normalizedImage, buffer);
            res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
        });
        svr.Get("/vis-terrain", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received GET request on /vis-terrain" << std::endl;
            std::vector<uchar> buffer;
            cv::imencode(".png", state->visualisation->terrainSavedImage, buffer);
            res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
        });

        svr.Post("/reset-simulation", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "Received POST request on /reset-simulation" << std::endl;
            state->simulation->reset();
            state->visualisation->simulationOffset = 0.0f;
            res.set_content("Simulation reset", "text/plain");
        });
        svr.Post("/run-simulation", [this](const httplib::Request&, httplib::Response& res) {
            std::lock_guard<std::mutex> lock(serverMutex);
            std::cout << "--------------------------" << std::endl;
            std::cout << "Running Simulation" << std::endl;
            state->visualisation->paused = true;
            state->simulation->saveToNetCDF(state->visualisation->terrainImage);
            state->simulation->saveToNetCDF(state->waterCanvas->waterDepthMat, "water.nc");
            state->simulation->testLoadSequence();
            state->visualisation->simulationOffset = state->simulation->pastOffset;
            res.set_content("Simulation run", "text/plain");
        });

        std::cout << "-------------------------" << std::endl;
        std::cout << "Starting server on 127.0.0.1:18080" << std::endl;
        svr.listen("localhost", 18080);
        std::cout << "Server stopped" << std::endl;
    }
    Remote(WindowState *windowState) : serverThread([this]() { startServer(); }), state(windowState) {
        serverThread.detach();
    }
};

#endif // REMOTE_HPP