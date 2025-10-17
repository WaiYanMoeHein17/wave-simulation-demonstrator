#include <glad/glad.h>
#include <opencv2/opencv.hpp>

class Cleaner {
public:
    GLuint texture = -1;

    cv::Mat kernel;
    cv::Mat closed;

    cv::Mat near_black_mask;
    int nearBlackThreshold = 150;

    cv::Mat dilated_mask;
    cv::Mat inpainted_image;
    cv::Mat grayscale_image;
    cv::Mat blurred_grayscale_image;

    // Cleaner() {};
    // ~Cleaner() {};

    void clean(const cv::Mat &cropped_image)
    {
        kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(12, 12));

        cv::morphologyEx(cropped_image, closed, cv::MORPH_CLOSE, kernel);

        cv::inRange(closed, cv::Scalar(0, 0, 0), cv::Scalar(nearBlackThreshold, nearBlackThreshold, nearBlackThreshold), near_black_mask);
        
        cv::dilate(near_black_mask, dilated_mask, kernel);
        
        cv::inpaint(closed, dilated_mask, inpainted_image, 1, cv::INPAINT_TELEA);

        cv::cvtColor(inpainted_image, grayscale_image, cv::COLOR_BGR2GRAY);

        cv::GaussianBlur(grayscale_image, blurred_grayscale_image, cv::Size(5, 5), 0);
    }

    void demoWindow(const cv::Mat &cropped_image)
    {
        cv::Mat bgc_image;
        cv::cvtColor(blurred_grayscale_image, bgc_image, cv::COLOR_GRAY2BGR);

        cv::Mat combined_top;
        cv::hconcat(cropped_image, closed, combined_top);

        cv::Mat combined_bottom;
        cv::hconcat(inpainted_image, bgc_image, combined_bottom);

        cv::Mat combined;
        cv::vconcat(combined_top, combined_bottom, combined);

        cv::namedWindow("Combined Image", cv::WINDOW_AUTOSIZE);
        cv::imshow("Combined Image", combined);
        cv::waitKey(0);
    }

    void toGL()
    {
        if (texture == -1) {
            glGenTextures(1, &texture);
        }
        glBindTexture(GL_TEXTURE_2D, texture);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, flipped.cols, flipped.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, flipped.data);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, blurred_grayscale_image.cols, blurred_grayscale_image.rows, 0, GL_RED, GL_UNSIGNED_BYTE, blurred_grayscale_image.data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glBindTexture(GL_TEXTURE_2D, 0);
    }
};