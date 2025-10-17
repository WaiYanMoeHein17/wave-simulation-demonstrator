#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <librealsense2/rs.hpp>
#include <librealsense2/rs_advanced_mode.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>

class Camera {
private:
public:
    bool isBagFile = false;
    bool playing = false;

    rs2::pipeline pipe;
    rs2::config cfg;
    rs2::colorizer colorMap;
    rs2::pipeline_profile profile;
    std::unique_ptr<rs2::playback> playback;

    rs2::align alignToDepth = rs2::align(RS2_STREAM_DEPTH);
    rs2::align alignToColor = rs2::align(RS2_STREAM_COLOR);

    cv::Size colorSize;
    cv::Size depthSize;

    Camera(std::string bagFile) {
        std::cout << "Initializing camera...\n";
        if (!bagFile.empty()) {
            isBagFile = true;
            std::cout << "Playing from bag file: " << bagFile << "\n";
            cfg.enable_device_from_file(bagFile);
        } else {
            std::cout << "Using RealSense device\n";
            cfg.enable_stream(RS2_STREAM_COLOR, 640, 480, RS2_FORMAT_BGR8, 30);
            cfg.enable_stream(RS2_STREAM_DEPTH, 848, 480, RS2_FORMAT_Z16, 30);
            // cfg.enable_stream(RS2_STREAM_DEPTH, 640, 480, RS2_FORMAT_Z16, 30);
            // cfg.enable_stream(RS2_STREAM_DEPTH, 1280, 720, RS2_FORMAT_Z16, 30);
            std::cout << "Enabled color and depth streams.\n";
        }

        if (!isBagFile) {
            rs2::context ctx;
            if (ctx.query_devices().size() == 0) {
                throw std::runtime_error("No RealSense devices found.");
            }
        }

        profile = pipe.start(cfg);
        std::cout << "Camera started.\n";

        if(isBagFile) {
            playback = std::make_unique<rs2::playback>(profile.get_device().as<rs2::playback>());
            playback->pause();
            playback->seek(std::chrono::milliseconds(0));
        }

        auto colorStream = profile.get_stream(RS2_STREAM_COLOR);
        auto depthStream = profile.get_stream(RS2_STREAM_DEPTH);

        colorSize = cv::Size(colorStream.as<rs2::video_stream_profile>().width(),
                             colorStream.as<rs2::video_stream_profile>().height());
        depthSize = cv::Size(depthStream.as<rs2::video_stream_profile>().width(),
                             depthStream.as<rs2::video_stream_profile>().height());

        std::cout << "Color Size: " << colorSize.width << "x" << colorSize.height << "\n";
        std::cout << "Depth Size: " << depthSize.width << "x" << depthSize.height << "\n";
        std::cout << "Camera initialized.\n";
    }

    void resume() {
        if (isBagFile)  playback->resume();
    }

    void warmUp() {
        if(isBagFile) return;
        std::cout << "Warming up the camera...\n";
        for (int i = 0; i < 30; ++i) {
            rs2::frameset frames;
            if (pipe.poll_for_frames(&frames)) {
                rs2::frame depthFrame = frames.get_depth_frame();
                rs2::frame colorFrame = frames.get_color_frame();
            }
        }
        std::cout << "Camera warmed up.\n";
    }

    void checkStoppedRestart() {
        if(!isBagFile) return;
        if (playing && playback->current_status() == RS2_PLAYBACK_STATUS_STOPPED) {
            std::cout << "Playback finished. Restarting...\n";
            pipe.stop();
            pipe.start(cfg);
            playback = std::make_unique<rs2::playback>(profile.get_device().as<rs2::playback>());
            playing = false;
        }
    }

    void softReset() {
        if (isBagFile) {
            playback->pause();
            playback->seek(std::chrono::milliseconds(0));
        }
    }

    void setROI(cv::Rect inROI) {
        if (isBagFile) {
            std::cout << "Cannot set ROI in playback mode.\n";
            return;
        }
        rs2::device dev = pipe.get_active_profile().get_device();
        rs2::sensor sensor = dev.first<rs2::depth_sensor>();
        rs2::roi_sensor s(sensor);
        // rs2::region_of_interest roi = s.get_region_of_interest();
        // s.set_region_of_interest(roi);
        std::cout << "ROI set to: " << inROI.x << ", " << inROI.y << ", " << inROI.width << ", " << inROI.height << "\n";
    }
};

#endif // CAMERA_HPP