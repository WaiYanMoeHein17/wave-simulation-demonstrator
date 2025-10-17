#ifndef SCREEN_HPP
#define SCREEN_HPP

#include <glm/glm.hpp>
#include <opencv2/opencv.hpp>

#include "Shader.hpp"

const char* basicVertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;

out vec2 TexCoord;

void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
    TexCoord = aTexCoord;
}
)";

const char* basicFragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

uniform sampler2D textureMap;

void main() {
    FragColor = texture(textureMap, TexCoord);
}
)";

class Screen {
private:
public:
    cv::Mat image;
    GLuint textureID;
    bool visible;


    GLuint VAO, VBO, shaderProgram;
    float quadVertices[24] = {
        // Positions     // Texture Coords
        -1.0f,  1.0f,   0.0f, 1.0f,  // Top-left
        -1.0f, -1.0f,   0.0f, 0.0f,  // Bottom-left
         1.0f, -1.0f,   1.0f, 0.0f,  // Bottom-right
    
        -1.0f,  1.0f,   0.0f, 1.0f,  // Top-left
         1.0f, -1.0f,   1.0f, 0.0f,  // Bottom-right
         1.0f,  1.0f,   1.0f, 1.0f   // Top-right
    };

    void initGL() {
        // Compile shaders
        unsigned int vertexShader = compileShader(GL_VERTEX_SHADER, basicVertexShaderSource);
        unsigned int fragmentShader = compileShader(GL_FRAGMENT_SHADER, basicFragmentShaderSource);
        
        // Create shader program
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);


        // Flip quadVertices vertically
        for (int i = 1; i < 24; i += 4) {
            quadVertices[i] = -quadVertices[i];
        }

        // Create VAO and VBO
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        
        glBindVertexArray(VAO);
        
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
        
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    Screen() : textureID(0), VAO(0), VBO(0), shaderProgram(0), visible(true) {
        initGL();
    }
    Screen(const cv::Mat& img) : textureID(0), VAO(0), VBO(0), shaderProgram(0), visible(true) {
        initGL();
        update(img);
    }
    ~Screen() {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
    }

    void draw() {
        if (!visible) return;
        glUseProgram(shaderProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glUniform1i(glGetUniformLocation(shaderProgram, "textureMap"), 0);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }

    void update(const cv::Mat& img) {
        image = img;


        if (textureID == 0)  glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);

        // texture params
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // to gl
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        if (img.channels() == 3) {
            // RGB image
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.cols, img.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, img.data);
        } else if (img.channels() == 1) {
            // Grayscale image
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, img.cols, img.rows, 0, GL_RED, GL_UNSIGNED_BYTE, img.data);
        } else if (img.depth() == CV_16U) {
            // Z16 image (16-bit depth)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R16, img.cols, img.rows, 0, GL_RED, GL_UNSIGNED_SHORT, img.data);
        } else {
            std::cerr << "Unsupported terrain format: channels=" << img.channels() 
                  << ", depth=" << img.depth() << std::endl;
        }
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, blurred_grayscale_image.cols, blurred_grayscale_image.rows, 0, GL_RED, GL_UNSIGNED_BYTE, blurred_grayscale_image.data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glBindTexture(GL_TEXTURE_2D, 0);
    }
};

#endif // VISUALISATION_HPP