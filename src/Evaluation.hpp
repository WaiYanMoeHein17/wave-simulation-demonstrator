#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include <opencv2/opencv.hpp>
#include <librealsense2/rs.hpp>

#include <string>
#include <vector>

double calculatePSNR(const cv::Mat& original, const cv::Mat& processed) {
    cv::Mat diff;
    cv::absdiff(original, processed, diff);
    diff.convertTo(diff, CV_32F);
    diff = diff.mul(diff);

    double mse = cv::sum(diff)[0] / (original.total() * original.channels());
    if (mse == 0) {
        return std::numeric_limits<double>::infinity(); // No difference
    }

    double psnr = 10.0 * std::log10((255 * 255) / mse);
    return psnr;
}

// https://docs.opencv.org/4.x/d5/dc4/tutorial_video_input_psnr_ssim.html
double calculateSSIM(const cv::Mat& original, const cv::Mat& processed) {
    const double C1 = 6.5025, C2 = 58.5225;

    cv::Mat originalFloat, processedFloat;
    original.convertTo(originalFloat, CV_32F);
    processed.convertTo(processedFloat, CV_32F);

    cv::Mat mu1, mu2;
    cv::GaussianBlur(originalFloat, mu1, cv::Size(11, 11), 1.5);
    cv::GaussianBlur(processedFloat, mu2, cv::Size(11, 11), 1.5);

    cv::Mat mu1Sq = mu1.mul(mu1);
    cv::Mat mu2Sq = mu2.mul(mu2);
    cv::Mat mu1Mu2 = mu1.mul(mu2);

    cv::Mat sigma1Sq, sigma2Sq, sigma12;
    cv::GaussianBlur(originalFloat.mul(originalFloat), sigma1Sq, cv::Size(11, 11), 1.5);
    sigma1Sq -= mu1Sq;

    cv::GaussianBlur(processedFloat.mul(processedFloat), sigma2Sq, cv::Size(11, 11), 1.5);
    sigma2Sq -= mu2Sq;

    cv::GaussianBlur(originalFloat.mul(processedFloat), sigma12, cv::Size(11, 11), 1.5);
    sigma12 -= mu1Mu2;

    cv::Mat t1 = 2 * mu1Mu2 + C1;
    cv::Mat t2 = 2 * sigma12 + C2;
    cv::Mat t3 = t1.mul(t2);

    t1 = mu1Sq + mu2Sq + C1;
    t2 = sigma1Sq + sigma2Sq + C2;
    t1 = t1.mul(t2);

    cv::Mat ssimMap;
    cv::divide(t3, t1, ssimMap);
    return cv::mean(ssimMap)[0];
}

class Evaluation {
private:
    cv::Mat groundTruth;
    cv::Mat image;
public:
    double accMean = 0.0;
    double accStd = 0.0;
    double accPSNR = 0.0;
    double accCount = 0.0;
    cv::Rect roi;

    void reset() {
        accMean = 0.0;
        accStd = 0.0;
        accPSNR = 0.0;
        accCount = 0.0;
        std::cout << "Accumulated metrics reset." << std::endl;
    }

    void accumulate(const cv::Scalar& mean, const cv::Scalar& stddev, const cv::Mat& img, const cv::Mat& prev) {
        accMean += mean[0];
        accStd += stddev[0];
        accPSNR += calculatePSNR(img, prev);
        accCount += 1.0;
        std::cout << "Mean: " << accMean / accCount << ", StdDev: " << accStd / accCount << ", PSNR: " << accPSNR / accCount << std::endl;
    }

    void setGroundTruth(const rs2::frame& depthFrame) {
        auto depthFrameData = depthFrame.get_data();
        auto width = depthFrame.as<rs2::video_frame>().get_width();
        auto height = depthFrame.as<rs2::video_frame>().get_height();

        groundTruth = cv::Mat(cv::Size(width, height), CV_16UC1, (void*)depthFrameData, cv::Mat::AUTO_STEP);
        groundTruth = groundTruth(roi);
    }
    void setImage(const rs2::frame& depthFrame) {
        auto depthFrameData = depthFrame.get_data();
        auto width = depthFrame.as<rs2::video_frame>().get_width();
        auto height = depthFrame.as<rs2::video_frame>().get_height();

        image = cv::Mat(cv::Size(width, height), CV_16UC1, (void*)depthFrameData, cv::Mat::AUTO_STEP);
        image = image(roi);
    }
    void evaluateCleaning() {
        if (image.empty()) {
            std::cerr << "Image or ground truth is empty." << std::endl;
            return;
        }
        cv::Mat normalizedImage;
        cv::normalize(image, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_8UC1);
        cv::Mat jetImage;
        cv::applyColorMap(normalizedImage, jetImage, cv::COLORMAP_JET);


        // gaussian blur
        cv::Mat gaussianBlurred;
        cv::GaussianBlur(image, gaussianBlurred, cv::Size(5, 5), 0);
        cv::Mat normalizedGaussian;
        cv::normalize(gaussianBlurred, normalizedGaussian, 0, 255, cv::NORM_MINMAX, CV_8UC1);
        cv::Mat jetGaussian;
        cv::applyColorMap(normalizedGaussian, jetGaussian, cv::COLORMAP_JET);

        // inpainting
        cv::Mat mask = (image == 0);
        cv::Mat inpainted;
        cv::inpaint(image, mask, inpainted, 3, cv::INPAINT_TELEA);
        cv::Mat normalizedInpainted;
        cv::normalize(inpainted, normalizedInpainted, 0, 255, cv::NORM_MINMAX, CV_8UC1);
        cv::Mat jetInpainted;
        cv::applyColorMap(normalizedInpainted, jetInpainted, cv::COLORMAP_JET);

        // dilation
        cv::Mat dilated;
        cv::dilate(image, dilated, cv::Mat(), cv::Point(-1, -1), 2);
        cv::Mat normalizedDilated;
        cv::normalize(dilated, normalizedDilated, 0, 255, cv::NORM_MINMAX, CV_8UC1);
        cv::Mat jetDilated;
        cv::applyColorMap(normalizedDilated, jetDilated, cv::COLORMAP_JET);
        
        // stack for display
        cv::Mat stackedImages;
        std::vector<cv::Mat> images = {jetImage, jetGaussian, jetInpainted, jetDilated};
        cv::hconcat(images, stackedImages);
        cv::imshow("Jet Images", stackedImages);

        groundTruth = image.clone();
        cv::Mat normalizedGroundTruth;
        cv::normalize(groundTruth, normalizedGroundTruth, 0, 255, cv::NORM_MINMAX, CV_8UC1);

        // calculate SSIM between ground truth and processed images
        double ssimOriginal = calculateSSIM(normalizedGroundTruth, normalizedImage);
        double ssimGaussian = calculateSSIM(normalizedGroundTruth, normalizedGaussian);
        double ssimInpainted = calculateSSIM(normalizedGroundTruth, normalizedInpainted);
        double ssimDilated = calculateSSIM(normalizedGroundTruth, normalizedDilated);
        // & print
        std::cout << "SSIM (Original): " << ssimOriginal << std::endl;
        std::cout << "SSIM (Gaussian): " << ssimGaussian << std::endl;
        std::cout << "SSIM (Inpainted): " << ssimInpainted << std::endl;
        std::cout << "SSIM (Dilated): " << ssimDilated << std::endl;


        cv::waitKey(0);
        cv::destroyAllWindows();
    }
};

#endif // EVALUATION_HPP