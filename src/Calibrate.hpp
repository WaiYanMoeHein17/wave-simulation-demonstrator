#ifndef CALIBRATE_HPP
#define CALIBRATE_HPP

#include <iostream>
#include <thread>

#include <librealsense2/rs.hpp>
#include <opencv2/opencv.hpp>

#include "Checkerboard.hpp"

class Calibrate {
private:
public:
    bool calibrated = false;
    std::vector<cv::Point2f> roiCorners;

    rs2::frameset getFirstFrames(rs2::pipeline& pipe) {
        rs2::frameset frames;
        while (!pipe.poll_for_frames(&frames)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        return frames;
    }

    void calibrate(rs2::pipeline& pipe) {
        // get first color & depth frames
        std::cout << "Calibrating..." << std::endl;
        auto frames = getFirstFrames(pipe);
        rs2::video_frame color_frame = frames.get_color_frame();
        rs2::depth_frame depth_frame = frames.get_depth_frame();
        
        // convert color frame to OpenCV format
        cv::Mat color_image(cv::Size(color_frame.get_width(), color_frame.get_height()), 
                    CV_8UC3, 
                    (void*)color_frame.get_data(), 
                    cv::Mat::AUTO_STEP);


        // Align color image to depth frame
        rs2::align align_to_depth(RS2_STREAM_DEPTH);
        frames = align_to_depth.process(frames);
        color_frame = frames.get_color_frame();
        cv::Mat aligned_color_image(cv::Size(color_frame.get_width(), color_frame.get_height()),
                              CV_8UC3,
                              (void*)color_frame.get_data(),
                              cv::Mat::AUTO_STEP);
        cv::cvtColor(aligned_color_image, aligned_color_image, cv::COLOR_RGB2BGR);
        color_image = aligned_color_image;
        // cv::cvtColor(color_image, color_image, cv::COLOR_RGB2BGR);

        // Select ROI & crop
        cv::Rect r = cv::selectROI("Select Sand Region", color_image);
        // cv::Rect r(
        //     (color_image.cols - 100) / 2,
        //     (color_image.rows - 100) / 2,
        //     100,
        //     100
        // );
        cv::destroyAllWindows();
        for (int i = 0; i < 5; ++i) {
            cv::waitKey(1);
        }
        cv::Mat roi = color_image(r);

        // ROI -> HSV
        cv::Mat hsv_roi;
        cv::cvtColor(roi, hsv_roi, cv::COLOR_BGR2HSV);

        // Calculate histogram for the ROI
        int histSize[] = { 30, 32 }; // Hue & Saturation
        // int histSize[] = { 60, 64 };
        
        float hranges[] = { 0, 180 };
        float sranges[] = { 0, 256 };
        const float* ranges[] = { hranges, sranges };
        int channels[] = { 0, 1 }; // Hue & Saturation
        cv::Mat hist;
        cv::calcHist(&hsv_roi, 1, channels, cv::Mat(), hist, 2, histSize, ranges, true, false);
        cv::normalize(hist, hist, 0, 255, cv::NORM_MINMAX);

        // full image -> HSV
        cv::Mat hsv_image;
        cv::cvtColor(color_image, hsv_image, cv::COLOR_BGR2HSV);

        // backprojection
        cv::Mat backproj;
        cv::calcBackProject(&hsv_image, 1, channels, hist, backproj, ranges, 10, true);
        cv::GaussianBlur(backproj, backproj, cv::Size(5, 5), 0);
        
        // threshold backprojection for mask
        cv::Mat mask;
        cv::threshold(backproj, mask, 25, 255, cv::THRESH_BINARY);

        // dilate
        const cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));
        // cv::dilate(mask, mask, kernel, cv::Point(-1, -1), 3);

        // cv::Mat blended;
        // cv::addWeighted(color_image, 0.7, overlay, 0.3, 0, blended);

        // find contours in the binary mask
        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(mask, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

        // largest contour
        auto largest_contour = *std::max_element(contours.begin(), contours.end(),
            [](const std::vector<cv::Point>& a, const std::vector<cv::Point>& b) {
                return cv::contourArea(a) < cv::contourArea(b);
            });

        // new mask w/ largest contour
        cv::Mat largest_contour_mask = cv::Mat::zeros(mask.size(), CV_8UC1);
        cv::drawContours(largest_contour_mask, std::vector<std::vector<cv::Point>>{largest_contour}, -1, cv::Scalar(255), cv::FILLED);

        // green overlay
        cv::Mat overlay = cv::Mat::zeros(color_image.size(), color_image.type());
        overlay.setTo(cv::Scalar(0, 255, 0), largest_contour_mask); // Set green color where the largest contour mask is non-zero

        
        if (!largest_contour.empty()) {
            // extreme points
            cv::Point top_left(mask.cols, mask.rows);
            cv::Point top_right(0, mask.rows);
            cv::Point bottom_left(mask.cols, 0);
            cv::Point bottom_right(0, 0);

            for (const auto& point : largest_contour) {
                if (point.x + point.y < top_left.x + top_left.y) top_left = point;
                if (point.x - point.y > top_right.x - top_right.y) top_right = point;
                if (point.x + point.y > bottom_right.x + bottom_right.y) bottom_right = point;
                if (point.x - point.y < bottom_left.x - bottom_left.y) bottom_left = point;
            }

            // overlay
            cv::addWeighted(color_image, 0.7, overlay, 0.3, 0, color_image);
            
            // draw extreme points
            cv::circle(color_image, top_left, 5, cv::Scalar(0, 0, 255), -1); // Red for top-left
            cv::circle(color_image, top_right, 5, cv::Scalar(0, 0, 255), -1); // Red for top-right
            cv::circle(color_image, bottom_left, 5, cv::Scalar(0, 0, 255), -1); // Red for bottom-left
            cv::circle(color_image, bottom_right, 5, cv::Scalar(0, 0, 255), -1); // Red for bottom-right

            roiCorners = { top_left, top_right, bottom_right, bottom_left };

            // debug
            cv::imshow("Rotated Rectangle Around Green Area", color_image);


            // save to timestampped file
            std::time_t now = std::time(nullptr);
            std::tm* local_time = std::localtime(&now);
            char filename[100];
            std::strftime(filename, sizeof(filename), "calibrated_image_%Y%m%d_%H%M%S.jpg", local_time);
            cv::imwrite(filename, color_image);
            cv::waitKey(0);
            cv::destroyAllWindows();
        } else {
            std::cerr << "No contours found!" << std::endl;
        }
    }

    void testDetectCheckerboard(Checkerboard& checkerboard, cv::Mat image) {
        // Load the image
        // cv::Mat image = cv::imread("/Users/macauley/Development/l4project/images/color_image_20250428_091050.png");
        // if (image.empty()) {
        //     std::cerr << "Failed to load image!" << std::endl;
        //     return;
        // }
        // image = cv::imread("/Users/macauley/Development/l4project/images/color_image_20250428_091050.png");

        // image -> gray
        cv::Mat gray_image;
        cv::cvtColor(image, gray_image, cv::COLOR_BGR2GRAY);

        // pattern size
        cv::Size pattern_size(7, 7);

        // corners?
        std::vector<cv::Point2f> corners;
        bool found = cv::findChessboardCorners(gray_image, pattern_size, corners,
                               cv::CALIB_CB_ADAPTIVE_THRESH | cv::CALIB_CB_NORMALIZE_IMAGE);

        if (found) {
            // refine locs
            cv::cornerSubPix(gray_image, corners, cv::Size(11, 11), cv::Size(-1, -1),
                     cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 30, 0.1));

            // debug draw on iamge
            cv::drawChessboardCorners(image, pattern_size, corners, found);

            // display for debug
            cv::imshow("Checkerboard Detection", image);
            cv::waitKey(0);
            cv::destroyAllWindows();

            // Compute the homography matrix
            if (corners.size() == checkerboard.points.size()) {
                cv::Mat homography = cv::findHomography(corners, checkerboard.points);

                std::cout << "Homography matrix:" << std::endl;
                std::cout << homography << std::endl;

                // warp image to fit & w/ perspective
                std::cout << "Warping image..." << std::endl;
                cv::Mat warped_image;
                cv::warpPerspective(image, warped_image, homography, checkerboard.image.size());



                // warp the extreme points
                std::cout << "Warping points:" << std::endl;
                std::vector<cv::Point2f> warped_corners;
                std::vector<cv::Point2f> roiCornersF = { roiCorners[0], roiCorners[1], roiCorners[2], roiCorners[3] };
                cv::perspectiveTransform(roiCornersF, warped_corners, homography);

                // draw
                for (const auto& point : warped_corners) {
                    std::cout << "Warped Point: (" << point.x << ", " << point.y << ")" << std::endl;
                    cv::circle(warped_image, point, 5, cv::Scalar(0, 0, 255), -1); // Red circles
                }

                // debug display
                cv::imshow("Warped Image", warped_image);
                cv::waitKey(0);
                cv::destroyAllWindows();
                checkerboard.update(warped_image);
            } else {
                std::cerr << "Mismatch between object points and detected corners!" << std::endl;
            }
        } else {
            std::cerr << "Checkerboard not found!" << std::endl;
        }

    }
};

#endif // CALIBRATE_HPP