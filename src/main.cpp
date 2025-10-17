/* 
L4 Advanced Project

Code Resources & Sources:
    *   https://learnopengl.com
        -   https://learnopengl.com/Getting-started/Creating-a-window
        -   https://learnopengl.com/Getting-started/Shaders
        -   https://learnopengl.com/Getting-started/Textures
        -   https://learnopengl.com/In-Practice/Text-Rendering
    *   https://docs.opencv.org/4.x/d5/dc4/tutorial_video_input_psnr_ssim.html
    *   https://github.com/IntelRealSense/librealsense
        -   https://github.com/IntelRealSense/librealsense/blob/master/doc/post-processing-filters.md
    *   https://web.cs.ucdavis.edu/~okreylos/ResDev/SARndbox/
    *   https://um-bridge-benchmarks.readthedocs.io/en/docs/
    *   https://github.com/annereinarz/ExaHyPE-Tsunami/tree/main/ApplicationExamples
    *   https://opencv.org
    *   https://www.unidata.ucar.edu/software/netcdf/examples/programs/
        -   https://www.unidata.ucar.edu/software/netcdf/examples/programs/SimpleXyWr.cpp
    *   https://examples.vtk.org/site/Cxx/
    *   https://threejs.org/examples/webgl_geometry_terrain_raycast.html
    *   ParaView
    *   GitHub Copilot (GPT-4o) vscode extension used for trivial code assist, auto-completion and snippets
*/


#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <thread>

#include <librealsense2/rs.hpp>
#include <opencv2/opencv.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <umbridge/umbridge.h>

#include <httplib.h>

#include "Visualisation.hpp"
#include "Window.hpp"
#include "Simulation.hpp"
#include "Calibrate.hpp"
#include "TextRenderer.hpp"
#include "Water.hpp"
#include "Checkerboard.hpp"
#include "Camera.hpp"
#include "Remote.hpp"
#include "Evaluation.hpp"


// #define ENABLE_EVALUATION 1

std::vector<cv::Point> selectPoints(const cv::Mat& image) {
    std::vector<cv::Point> points;
    cv::namedWindow("Select Points", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Select Points", [](int event, int x, int y, int, void* userdata) {
        if (event == cv::EVENT_LBUTTONDOWN) {
            std::vector<cv::Point>* points = reinterpret_cast<std::vector<cv::Point>*>(userdata);
            points->emplace_back(x, y);
            std::cout << "Point selected: (" << x << ", " << y << ")" << std::endl;
        }
    }, &points);

    while (points.size() < 4) {
        cv::imshow("Select Points", image);
        cv::waitKey(1);
    }

    cv::destroyWindow("Select Points");
    return points;
}


int main(int argc, char **argv) try
{
    // arguments
    std::string imageInputPath;
    std::string bagFile;
    std::vector<cv::Point> points;
    bool fullscreen = false;
    bool shouldCalibrate = false;
    std::string host = "http://localhost:4242";
    // std::string simulationInputPath = "/Users/macauley/Development/T/output2";
    std::string simulationInputPath = "/Users/macauley/Development/T/out";
    std::string simulationOutputPath = "/Users/macauley/Development/T/in";

    std::string diffPath;

    float temporalAlpha = 0.1f;
    float temporalDelta = 60.0f;

    // parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--image") {
            if (i + 1 < argc) {
            imageInputPath = argv[++i];
            } else {
                throw std::invalid_argument("No image path specified after --image");
            }
        }else if (arg == "--bag") {
            if (i + 1 < argc) {
                bagFile = argv[++i];
            } else {
                throw std::invalid_argument("No bag file specified after --bag");
            }
        } else if (arg == "--fullscreen") {
            fullscreen = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [--bag <file>] [--help]\n";
            exit(0);
        } else if (arg == "--points") {
            if (i + 4 < argc) {
                for (int j = 0; j < 4; ++j) {
                    int x, y;
                    std::stringstream ss(argv[++i]);
                    char comma;
                    if (!(ss >> x >> comma >> y) || comma != ',') {
                        throw std::invalid_argument("Invalid point format. Use x,y");
                    }
                    points.emplace_back(x, y);
                }
            } else {
            throw std::invalid_argument("Insufficient points specified after --points");
            }
        } else if (arg == "--host") {
            if (i + 1 < argc) {
            host = argv[++i];
            } else {
            throw std::invalid_argument("No host specified after --host");
            }
        } else if (arg == "--simulationInput") {
            if (i + 1 < argc) {
                simulationInputPath = argv[++i];
            } else {
                throw std::invalid_argument("No simulation input path specified after --simulationInput");
            }
        } else if (arg == "--simulationOutput") {
            if (i + 1 < argc) {
                simulationOutputPath = argv[++i];
            } else {
                throw std::invalid_argument("No simulation output path specified after --simulationOutput");
            }
        } else if (arg == "--calibrate") {
            shouldCalibrate = true;
        } else if (arg == "--diff") {
            if (i + 1 < argc) {
                diffPath = argv[++i];
            } else {
                throw std::invalid_argument("No diff path specified after --diffPath");
            }
        } else if (arg == "--temporalAlpha") {
            if (i + 1 < argc) {
            temporalAlpha = std::stof(argv[++i]);
            } else {
            throw std::invalid_argument("No value specified after --temporalAlpha");
            }
        } else if (arg == "--temporalDelta") {
            if (i + 1 < argc) {
            temporalDelta = std::stof(argv[++i]);
            } else {
            throw std::invalid_argument("No value specified after --temporalDelta");
            }
        } else {
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }
    bool imageMode = !imageInputPath.empty();

    // https://github.com/UM-Bridge/umbridge/blob/main/clients/c%2B%2B/http-client.cpp
    // std::cout << "Connecting to host " << host << std::endl;

    // // List supported models
    // std::vector<std::string> models = umbridge::SupportedModels(host);
    // std::cout << "Supported models: " << std::endl;
    // for (auto model : models) {
    // std::cout << "  " << model << std::endl;
    // }
    // umbridge::HTTPModel client(host, "forward");
    // std::vector<std::vector<double>> inputs {{100.0, 18.0}};
    // std::vector<std::vector<double>> outputs = client.Evaluate(inputs);
    // std::cout << "Output: ";
    // for (const auto& value : outputs[0]) {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;

    Simulation sim(simulationInputPath, simulationOutputPath, host);

    Camera camera(bagFile);

    cv::Mat image;
    if (imageMode) {
        std::cout << "Using image input: " << imageInputPath << "\n";
        cv::Mat image = cv::imread(imageInputPath);
        if (image.empty()) {
            std::cerr << "Failed to load image: " << imageInputPath << "\n";
            return -1;
        }
        points = selectPoints(image);

        // crop to selection
        cv::Rect boundingBox = cv::boundingRect(points);
        cv::Mat croppedDepthMat = image(boundingBox);

        // jet -> gray
        cv::Mat redChannel(croppedDepthMat.size(), CV_8UC1);
        int fromTo[] = { 2, 0 };
        cv::mixChannels(&croppedDepthMat, 1, &redChannel, 1, fromTo, 1);

        cv::Mat grayscaleDepthMat = redChannel.clone();

        cv::imshow("Cropped Depth Map", grayscaleDepthMat);
        cv::waitKey(0);

        // image = invertedDepthMat;
    }
    Calibrate calibrate;
    if (shouldCalibrate) {
        std::cout << "-------------------------------------------\n";
        std::cout << "Calibrating camera...\n";
        camera.resume();
        calibrate.calibrate(camera.pipe);
        camera.softReset();
        cv::Point2f center(0, 0);
        for (const auto& point : calibrate.roiCorners) {
            center += cv::Point2f(point);
        }
        center *= (1.0f / calibrate.roiCorners.size());

        for (const auto& point : calibrate.roiCorners) {
            cv::Point2f vec = cv::Point2f(point) - center;
            cv::Point2f scaled = center + vec * 0.90f; // scale in by 5%
            points.push_back(cv::Point(cvRound(scaled.x), cvRound(scaled.y)));
        }
        std::cout << "ROI corners: " << points[0] << ", " << points[1] << ", " << points[2] << ", " << points[3] << "\n";

        std::cout << "Settting Camera ROI" << std::endl;
        cv::Rect roiRect = cv::boundingRect(points);
        // camera.setROI(roiRect);
    }

    // cv::Rect rrr(0, 0, camera.colorSize.width, camera.colorSize.height);
    // camera.setROI(rrr);


    std::cout << "-------------------------------------------\n";
    std::cout << "Window Init...\n";
    // GLFW Window INIT
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    // https://discourse.glfw.org/t/fullscreen-minimizes-on-clicking-on-another-monitor/1451
    glfwWindowHint( GLFW_AUTO_ICONIFY, GLFW_FALSE );

    int count;
    GLFWmonitor** monitors = glfwGetMonitors(&count);
    std::cout << "Number of monitors: " << count << std::endl;
    GLFWmonitor* monitor = count > 1 ? monitors[1] : monitors[0];
    const GLFWvidmode* mode = glfwGetVideoMode(monitor);
    if (!mode) {
        std::cerr << "Failed to get video mode\n";
        glfwTerminate();
        return -1;
    }
    std::cout << "Monitor resolution: " << mode->width << "x" << mode->height << std::endl;

    // Create a GLFW window
    // GLFWwindow* window = glfwCreateWindow(800, 600, "GL", nullptr, nullptr);
    GLFWwindow* window = fullscreen ? glfwCreateWindow(mode->width, mode->height, "GL", monitor, nullptr) : glfwCreateWindow(800, 600, "GL", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Load OpenGL functions using GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  

    // Enable anti-aliasing
    // glfwWindowHint(GLFW_SAMPLES, 4);
    // glEnable(GL_MULTISAMPLE);

    // Set up viewport
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    std::cout << "Frame Buffer size: " << width << "x" << height << std::endl;

    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);

    TextRenderer textRenderer("../shaders/text.vs", "../shaders/text.fs", windowWidth, windowHeight);

    Visualisation vis;
    Water water(windowWidth, windowHeight);

    WindowState windowState;
    windowState.visualisation = &vis;
    windowState.simulation = &sim;
    windowState.waterCanvas = &water;
    glfwSetWindowUserPointer(window, &windowState);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    
    //-------------------------------------
    // cursor
    cv::Mat cursor(50, 50, CV_8UC4, cv::Scalar(0, 0, 0, 0));
    int thickness = 4;
    cv::line(cursor, cv::Point(25, 10), cv::Point(25, 40), cv::Scalar(255, 255, 255, 255), thickness); // y
    cv::line(cursor, cv::Point(10, 25), cv::Point(40, 25), cv::Scalar(255, 255, 255, 255), thickness); // x
    GLFWimage cursorImage;
    cursorImage.width = cursor.cols;
    cursorImage.height = cursor.rows;
    cursorImage.pixels = cursor.data;
    GLFWcursor* customCursor = glfwCreateCursor(&cursorImage, 25, 25); // center
    if (customCursor) {
        glfwSetCursor(window, customCursor);
    } else {
        std::cerr << "Failed to create custom cursor" << std::endl;
    }
    //-------------------------------------

    if(!diffPath.empty()) {
        std::cout << "Diff path: " << diffPath << std::endl;
        windowState.differenceCanvas = new Difference(diffPath);
    }

    Remote remote(&windowState);

    // Checkerboard checkerboard(1280, 720);
    Checkerboard checkerboard(fullscreen ? mode->width : windowWidth, 
                               fullscreen ? mode->height : windowHeight);
    
    camera.resume();
    // if (shouldCalibrate) {
    //     // checkerboard.debugCheckerboard();
    //     rs2::frameset frameset;
    //     if (camera.pipe.poll_for_frames(&frameset)) {
    //         rs2::frame colorFrame = frameset.get_color_frame();
    //         cv::Mat colorMat(cv::Size(colorFrame.as<rs2::video_frame>().get_width(),
    //                                   colorFrame.as<rs2::video_frame>().get_height()),
    //                          CV_8UC3, (void*)colorFrame.get_data(), cv::Mat::AUTO_STEP);
    //         cv::cvtColor(colorMat, colorMat, cv::COLOR_RGB2BGR);
    //         // cv::imshow("Color Frame", colorMat);
    //         // cv::waitKey(1);
    //         // calibrate.testDetectCheckerboard(checkerboard, colorMat);
    //     }
    // }


    Evaluation eval;

    // Main loop
    double lastTime = glfwGetTime();
    double deltaTime = 0.0;
    
    cv::Mat prev;
    cv::Rect boundingBox = points.size() > 0 ? cv::boundingRect(points) : cv::Rect(0, 0, 0, 0);
    eval.roi = boundingBox;
    rs2::frameset frames;
    rs2::spatial_filter spatial_filter;
    rs2::temporal_filter temporal_filter;
    rs2::hole_filling_filter hole_filling;

    temporal_filter.set_option(RS2_OPTION_FILTER_SMOOTH_ALPHA, temporalAlpha); // higher alpha -> stronger smoothing
    temporal_filter.set_option(RS2_OPTION_FILTER_SMOOTH_DELTA, temporalDelta); // larger delta -> more aggressive filtering

    hole_filling.set_option(RS2_OPTION_HOLES_FILL, 2);

    double lastMotionTime = 0.0;
    // playback.resume();
    camera.warmUp();
    std::cout << "Before loop" << std::endl;
    while (!glfwWindowShouldClose(window)) {
        double currentTime = glfwGetTime();
        deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        camera.checkStoppedRestart();
        
        bool motionDetected = false;
        if (camera.pipe.poll_for_frames(&frames)) {
            // auto start = std::chrono::high_resolution_clock::now();
            camera.playing = true;

            rs2::frame colorFrame = frames.get_color_frame();
            rs2::frame depthFrame = frames.get_depth_frame();
            rs2::frame depthColorized = camera.colorMap.colorize(depthFrame);

            if (glfwGetKey(window, GLFW_KEY_F9) == GLFW_PRESS) {
                std::cout << "Setting ground truth..." << std::endl;
                eval.setGroundTruth(depthFrame);
            }
            if (glfwGetKey(window, GLFW_KEY_F10) == GLFW_PRESS) {
                std::cout << "Setting image..." << std::endl;
                eval.setImage(depthFrame);
                eval.evaluateCleaning();
            }

            // filters
            if(windowState.enableFilter) {
                // depthFrame = spatial_filter.process(depthFrame);
                depthFrame = temporal_filter.process(depthFrame);
                depthFrame = hole_filling.process(depthFrame);
            }

            if (glfwGetKey(window, GLFW_KEY_F7) == GLFW_PRESS) {
                std::cout << "Setting image AFTER RS..." << std::endl;
                eval.setImage(depthFrame);
                eval.evaluateCleaning();
            }

            cv::Mat depthMat(camera.depthSize, 
                        CV_16UC1, 
                        (void*)depthFrame.get_data(), 
                        cv::Mat::AUTO_STEP);

            cv::Mat depthColorizedMat(camera.depthSize, 
                        CV_8UC3, 
                        (void*)depthColorized.get_data(), 
                        cv::Mat::AUTO_STEP);

            if (glfwGetKey(window, GLFW_KEY_F5) == GLFW_PRESS) eval.reset();
            if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS) {
                std::time_t now = std::time(nullptr);
                std::stringstream filename;
                filename << "color_image_" << std::put_time(std::localtime(&now), "%Y%m%d_%H%M%S") << ".png";
                cv::Mat colorMat(cv::Size(colorFrame.as<rs2::video_frame>().get_width(), 
                                          colorFrame.as<rs2::video_frame>().get_height()), 
                                 CV_8UC3, (void*)colorFrame.get_data(), cv::Mat::AUTO_STEP);
                cv::imwrite(filename.str(), colorMat);
                std::cout << "Saved color image to " << filename.str() << std::endl;
            }

            // cv::Mat noralizedDepthMat;
            // cv::normalize(depthMat, noralizedDepthMat, 0, 255, cv::NORM_MINMAX, CV_8UC1);
            // cv::bitwise_not(noralizedDepthMat, noralizedDepthMat);
            // cv::imshow("Depth Map", noralizedDepthMat);

            if (points.empty()) {
                rs2::frameset aligned_frames = camera.alignToDepth.process(frames);
                rs2::frame aligned_color_frame = aligned_frames.get_color_frame();
                cv::Mat colorMat(cv::Size(aligned_color_frame.as<rs2::video_frame>().get_width(),
                                          aligned_color_frame.as<rs2::video_frame>().get_height()),
                                 CV_8UC3, (void*)aligned_color_frame.get_data(), cv::Mat::AUTO_STEP);
                cv::cvtColor(colorMat, colorMat, cv::COLOR_RGB2BGR);

                points = selectPoints(colorMat);
                // points = selectPoints(depthColorizedMat);
                boundingBox = cv::boundingRect(points);
                eval.roi = boundingBox;
                std::cout << "Bounding box: " << boundingBox << std::endl;
            }

            //  = cv::boundingRect(points);
            cv::Mat croppedDepthMat = depthMat(boundingBox);
            cv::Mat normalizedDepthMat;
            cv::normalize(croppedDepthMat, normalizedDepthMat, 0, 255, cv::NORM_MINMAX, CV_8UC1);
            cv::bitwise_not(normalizedDepthMat, normalizedDepthMat);
            // cv::imshow("Depth Colorized", normalizedDepthMat);

            if (windowState.saveNext) {
                std::time_t now = std::time(nullptr);
                std::stringstream filename;
                filename << "../results/" << std::put_time(std::localtime(&now), "%Y%m%d_%H%M%S") << ".png";
                cv::imwrite(filename.str(), depthColorizedMat(boundingBox));
                std::cout << "Saved depth colorized image to " << filename.str() << std::endl;
                windowState.saveNext = false;
            }
            
            if (!prev.empty()) {
                cv::Mat diff;
                cv::absdiff(normalizedDepthMat, prev, diff);
                cv::Mat threshDiff;
                cv::threshold(diff, threshDiff, 30, 255, cv::THRESH_BINARY);
                
                #ifdef ENABLE_EVALUATION
                cv::Scalar mean, stddev;
                cv::meanStdDev(diff, mean, stddev);
                eval.accumulate(mean, stddev, normalizedDepthMat, prev);
                #endif
                
                double motion = cv::sum(threshDiff)[0] / 255.0;
                motionDetected = motion > 100;
                if (motion > 100) {
                    // std::cout << "Motion detected: " << motion << " pixels" << std::endl;
                    lastMotionTime = currentTime;
                }

                if (windowState.debugWindows) {
                    cv::Mat jet;
                    cv::applyColorMap(normalizedDepthMat, jet, cv::COLORMAP_JET);
                    if (windowState.debugWindows && motion > 100) 
                        cv::putText(jet, "MOTION", cv::Point(10, 30), cv::FONT_HERSHEY_SIMPLEX, 
                                    0.5, cv::Scalar(0, 255, 0), 1);
                    cv::imshow("Original Depth Map", jet);
                    cv::imshow("Motion Mask", threshDiff);
                    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS) {
                        std::time_t now = std::time(nullptr);
                        std::stringstream motionFilename, colorFilename;

                        motionFilename << "../results/motion_mask_" << std::put_time(std::localtime(&now), "%Y%m%d_%H%M%S") << ".png";
                        colorFilename << "../results/color_image_" << std::put_time(std::localtime(&now), "%Y%m%d_%H%M%S") << ".png";

                        cv::imwrite(motionFilename.str(), threshDiff);

                        cv::Mat colorMat(cv::Size(colorFrame.as<rs2::video_frame>().get_width(),
                                                  colorFrame.as<rs2::video_frame>().get_height()),
                                         CV_8UC3, (void*)colorFrame.get_data(), cv::Mat::AUTO_STEP);
                        cv::cvtColor(colorMat, colorMat, cv::COLOR_RGB2BGR);
                        cv::imwrite(colorFilename.str(), colorMat);

                        std::cout << "Saved motion mask to " << motionFilename.str() << std::endl;
                        std::cout << "Saved color image to " << colorFilename.str() << std::endl;
                    }
                }
            }
            prev = normalizedDepthMat.clone();

            if (currentTime - lastMotionTime > 1.0) {
                std::lock_guard<std::mutex> lock(remote.serverMutex);
        
                vis.terrainToGL(normalizedDepthMat);
            }
        }
        std::lock_guard<std::mutex> lock(remote.serverMutex);


        // update simulation texture every 1/2 second
        // static double lastTextureUpdateTime = 0.0;
        // if (currentTime - lastTextureUpdateTime >= 0.1) {
        //     sim.advanceFrame();
        //     lastTextureUpdateTime = currentTime;
        // }
        float simMaxValue = 0.0;
        sim.advanceFrame(simMaxValue);

        // if (water.maxValue > 1.0f) {
        //     vis.simulationScale = water.maxValue / 2.55f;
        //     std::cout << "Max value: " << water.maxValue << std::endl;
        //     std::cout << "Simulation scale: " << vis.simulationScale << std::endl;
        // }
        // vis.draw(sim.texture, water.texture);


        bool simOrWater = sim.frameCount() == 0;
        vis.simulationScale = simOrWater ? 2.55 : simMaxValue;
        vis.draw((simOrWater) ? water.texture : sim.texture, 
                 (simOrWater) ? &water.waterTextureMat : sim.getCurrentTexture(),
                 !simOrWater);

        
        if(windowState.markers.size() == 2 && !simOrWater) {
            vis.samplePointsAlongLine(windowState.markers[0], windowState.markers[1], 100, 
                windowWidth, windowHeight,
                *sim.getCurrentFrame(),  sim.currentFrame - 1
            );
            
        }

        if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) {
            checkerboard.draw();
        }
        if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS && !calibrate.calibrated) {
            rs2::frame colorFrame = frames.get_color_frame();
            cv::Mat colorMat(cv::Size(colorFrame.as<rs2::video_frame>().get_width(), 
            colorFrame.as<rs2::video_frame>().get_height()), 
            CV_8UC3, (void*)colorFrame.get_data(), cv::Mat::AUTO_STEP);
            calibrate.testDetectCheckerboard(checkerboard, colorMat);
            calibrate.calibrated = true;
            glfwFocusWindow(window);
            std::cout << "Calibration complete?" << std::endl;
        }
        if(glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
            calibrate.calibrated = false;
            checkerboard.reset();
            std::cout << "Calibration reset" << std::endl;
        }
        if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS && windowState.differenceCanvas != nullptr) {
            windowState.differenceCanvas->draw();
        }



        std::stringstream statusText;
        if(camera.isBagFile) {
            auto playback_duration = std::chrono::duration_cast<std::chrono::seconds>(camera.playback->get_duration()).count();
            uint64_t current_time = camera.playback->get_position() / 1e9;    
            statusText << "Playback: " << current_time << " / " << playback_duration << " seconds";
        }
        // statusText << " | Sim. Frame: " << sim.currentFrame << " / " << sim.frames.size();
        if (vis.isPaused()) statusText << " | (Vis. Paused)";

        if (sim.isRunning) {
            statusText << " | (Sim. Running)";
        } else{
            statusText << " | Sim. Frame: " << sim.currentFrame << " / " << sim.frameCount();
        }

        if(motionDetected) {
            statusText << " | Motion Detected";
        }

        statusText << " | Sim. Scale: " << vis.simulationScale;

        textRenderer.renderText(statusText.str(), 10.0f, 10.0f, 0.25f, glm::vec3(1.0f, 1.0f, 1.0f));

        if (sim.isRunning) {
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            // static double simStartTime = glfwGetTime();
            double elapsed = glfwGetTime() - sim.simStartTime;
            double remaining = 13.0 - elapsed;
            if (remaining < 0.0) remaining = 0.0;

            float centerX = windowWidth / 2.0f;
            float centerY = windowHeight / 2.0f;
            
            std::string simText = "Simulation Running";
            std::stringstream timerText;
            timerText << "Est. Time left: " << std::fixed << std::setprecision(1) << remaining << " s";

            textRenderer.renderText(
                simText,
                centerX - 200.0f,
                centerY,
                1.0f,
                glm::vec3(0.0f, 0.0f, 0.0f)
            );
            textRenderer.renderText(
                timerText.str(),
                centerX - 150.0f,
                centerY + 60.0f,
                0.8f,
                glm::vec3(0.0f, 0.0f, 0.0f)
            );
        }else{
            sim.simStartTime = glfwGetTime();
        }

        // terrain image for api every 5 seconds
        static double lastSaveTime = 0.0;
        if (currentTime - lastSaveTime >= 5.0) {
            vis.terrainSavedImage = vis.getGLImage();
            lastSaveTime = currentTime;
        }

        // Swap buffers
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

    camera.pipe.stop();
    cv::destroyAllWindows();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
} catch (const rs2::error & e) {
    std::cerr << "RealSense error: " << e.what() << "\n";
    return EXIT_FAILURE;
} catch (const std::exception & e) {
    std::cerr << "Exception: " << e.what() << "\n";
    return EXIT_FAILURE;
}
