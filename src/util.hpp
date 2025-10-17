#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

void cvMat2GLtexture(const cv::Mat& image, GLuint& texture) {
    if (texture == 0) {
        glGenTextures(1, &texture);
    }
    glBindTexture(GL_TEXTURE_2D, texture);

    // texture params
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // image data to GL
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    if (image.channels() == 3) {
        // RGB image
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image.cols, image.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, image.data);
    } else if (image.channels() == 1) {
        // Grayscale image
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, image.cols, image.rows, 0, GL_RED, GL_UNSIGNED_BYTE, image.data);
    } else if (image.depth() == CV_16U) {
        // Z16 image (16-bit depth)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16, image.cols, image.rows, 0, GL_RED, GL_UNSIGNED_SHORT, image.data);
    } else {
        std::cerr << "Unsupported image format: channels=" << image.channels() 
              << ", depth=" << image.depth() << std::endl;
    }
    // glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, blurred_grayscale_image.cols, blurred_grayscale_image.rows, 0, GL_RED, GL_UNSIGNED_BYTE, blurred_grayscale_image.data);
    glGenerateMipmap(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);
}

#endif // UTIL_HPP